function F = compute_features_from_wave(f, init, beta, gamma, params)
% COMPUTE_FEATURES_FROM_WAVE (hybrid-threshold, hardened)
% - 列向量/清洗/短序列保护
% - t_cross: 采用分叉阈值（负相紧、近零松）+ 负相优先 + 连续m点（默认1）
% - 其它特征仍用 eps0 = c_eps * MADmean（保持你原有口径稳定）
% - κ：动态窗口 + 差分回退；尾部：主口径14天 + 相对口径

% ---- defaults ----
C = struct('c_eps',0.3,'q_early',0.30,'tail_ref',14, ...
           'c_eps_hi',0.3,'c_eps_lo',0.3, ...
           'alpha_hi',0.1,'alpha_lo',0.10, ...
           'eps_mode',"hybrid-min", ...  % 供将来扩展；当前分叉逻辑不依赖此开关
           'm_neg',1,'m_nz',1);
if nargin>=5 && ~isempty(params)
    fn = fieldnames(params);
    for i=1:numel(fn), C.(fn{i}) = params.(fn{i}); end
end
C.eps_mode = string(C.eps_mode);

% ---- 强制列向量 & 清洗 ----
f = double(f(:));
T = numel(f); t = (1:T)';
F = struct(); F.T = T;
if T==0, error('empty f provided'); end
bad = ~isfinite(f);
if any(bad)
    if all(bad), f(:)=0;
    else
        f = fillmissing(f,'linear','EndValues','nearest');
        f(~isfinite(f)) = 0;
    end
end

% ---- 基础 eps0（仅供除 t_cross 以外的特征使用；保持原口径）----
mu   = mean(f);
mad_mean = mean(abs(f - mu));
eps0 = C.c_eps * mad_mean;
F.eps_val  = eps0;      % 老字段（面积/计数等仍用它）
F.eps_base = mad_mean;

% ---- t_cross 专用的分叉阈值 ----
fmax = max(f);
eps_hi  = min(C.c_eps_hi * mad_mean, C.alpha_hi * fmax);   % 正侧上阈（紧）
eps_neg = min(C.c_eps_lo * mad_mean, C.alpha_lo * fmax);   % 负相阈值（紧）
eps_nz  = max(C.c_eps_lo * mad_mean, C.alpha_lo * fmax);   % 近零阈值（松）
F.eps_hi  = eps_hi;
F.eps_neg = eps_neg;
F.eps_nz  = eps_nz;

% ---- ODE integration ----
traj = integrate_seir_odeint_style(f, init, beta, gamma);
S = traj.S(:); I = traj.I(:); R = traj.R(:);
minT = min([T, numel(S), numel(I), numel(R)]);
if minT < T
    f = f(1:minT); t = (1:minT)'; T = minT;
    S = S(1:minT); I = I(1:minT); R = R(1:minT);
end
N = S + I + R; N(N<=0) = eps;
I_ode = traj.I_ode(:); I_ode = I_ode(1:T);

% ==== 峰值/峰时刻 ====
[F.I_ode_peak, F.t_I_ode_peak] = safe_max(I_ode);
[F.S_peak,     F.t_S_peak    ] = safe_max(S);
[F.I_peak,     F.t_I_peak    ] = safe_max(I);
[F.R_peak,     F.t_R_peak    ] = safe_max(R);

% ---- 基本量 ----
inc  = beta * (S .* I) ./ N;
dinc = [0; diff(inc)];

% ---- 主波上升段内的 Δinc 峰 ----
[~, idx_fmax] = max(f); F.t_fmax = idx_fmax;
hi = max(2, idx_fmax);
[~, rel] = max(dinc(2:hi));
idx_inc_slope_main = rel + 1;
F.t_inc_slope_max = idx_inc_slope_main;

% 审计：全时段 Δinc 峰
[~, idx_all] = max(dinc);
F.t_inc_slope_max_all = max(1, idx_all);

% ---- crossing：正侧 -> 负(优先)/近零（分叉阈值 + 连续m点，默认m=1）----
% ---- crossing（统一：t_cross = 正侧穿越 +epsilon 的时刻）----
t_pos = NaN;    % earliest time dropping from >eps_hi to <= +eps_nz (m_nz consecutive)
t_neg = NaN;    % earliest time dropping from >eps_hi to <  -eps_neg (m_neg consecutive)
cross_type = -1;

if T > idx_fmax
    % 1) near-zero：来自正侧，连续 mN 天 ≤ eps_nz
    mN = max(1, C.m_nz);        % <<<<<< 改：允许 m_nz = 0 时退化为 1
    for tt = (idx_fmax+1):T
        if f(tt-1) > eps_hi && tt + mN - 1 <= T ...
                && all( f(tt:tt+mN-1) <= eps_nz )
            t_pos = tt; 
            break;
        end
    end

    % 2) negative：来自正侧，连续 mG 天 < -eps_neg
    mG = max(1, C.m_neg);       % <<<<<< 改：允许 m_neg = 0 时退化为 1
    for tt = (idx_fmax+1):T
        if f(tt-1) > eps_hi && tt + mG - 1 <= T ...
                && all( f(tt:tt+mG-1) <  -eps_neg )
            t_neg = tt; 
            break;
        end
    end
end

% 3) 分类与赋值：类型按是否负相；t_cross 一律取“正穿越 +epsilon”的时刻
if ~isnan(t_neg)
    cross_type = 1;     % 负相
elseif ~isnan(t_pos)
    cross_type = 0;     % 近零
else
    cross_type = -1;    % 删失
end

F.t_cross    = t_pos;       % 统一为正的穿越 +epsilon 的时刻（与 cross_type 无关）
F.cross_type = cross_type;

% ★★★ 为后续仍引用的局部变量赋值
t_cross = F.t_cross;

% ---- delta_early ----
if ~isnan(F.t_cross), F.delta_early = F.t_cross - F.t_inc_slope_max; else, F.delta_early = NaN; end
F.delta_early_rel = F.delta_early / T;

F.is_shortwave    = double(T < 14);

% ---- κ（动态窗口 + 回退）----
if ~isnan(t_cross) && T>=2
    w_dyn = max(1, min(3, floor(0.1 * T)));
    lo = max(1, t_cross - w_dyn);
    hi = min(T, t_cross + w_dyn);
    if hi-lo+1 >= 2
        p = polyfit((lo:hi)', f(lo:hi), 1);
        F.kappa_cross = p(1);
    elseif t_cross>1 && t_cross<T
        F.kappa_cross = (f(t_cross+1) - f(t_cross-1))/2;
    elseif t_cross>1
        F.kappa_cross = f(t_cross) - f(t_cross-1);
    elseif t_cross<T
        F.kappa_cross = f(t_cross+1) - f(t_cross);
    else
        F.kappa_cross = NaN;
    end
else
    F.kappa_cross = NaN;
end

% ---- 面积 & 比例（仍用 eps0；dt=1 用 sum）----
fpos = f >  eps0; 
fneg = f < -eps0;
A_plus       = sum( f(fpos) );
A_minus      = sum( abs(f(fneg)) );
A_minus_eps  = sum( max(0, eps0 - f(~fpos)) );
F.A_plus = A_plus; F.A_minus = A_minus; F.A_minus_eps = A_minus_eps;
F.rho_eps = (A_minus + 0.3*A_minus_eps) / max(A_plus, eps);

% ---- 无量纲 ----
tf1 = (beta ./ N) .* f;
early_len = max(1, floor(C.q_early * T));
F.tf1_mean_early = mean(tf1(1:early_len));
ratio = f ./ max(inc, eps);
[F.tf2_peak, F.t_tf2_peak] = safe_max(ratio);

% ==== 结构 & 尾部（更稳健的多波判定）====
% 轻度平滑（仅用于形状判定，不改变其它计算）
w_smooth = max(3, 2*floor(0.05*T)+1);   % ~5% 波长，奇数窗
f_s = movmean(f, w_smooth, 'omitnan');

% 峰值筛选参数（内部阈值，足够保守避免误检）
minPeakHeight   = max(eps0, 0.3*max(f_s));    % 峰高至少超过近零并占主峰的10%
minPeakDist     = max(7, round(0.20*T));       % 峰间至少相隔~10%波长（且≥7天）
promRel         = 0.3;                        % 峰突出度至少占全局幅度的25%
valleyDropRel   = 0.2;                        % 两峰之间谷底至少下降到 min(两峰)*50% 或低于 eps0

% 取峰：优先 findpeaks；若不可用则退化到手写简单峰再做阈值过滤
pk = []; loc = []; prom = [];
if exist('findpeaks','file') == 2
    [pk, loc, w_pk, prom] = findpeaks(f_s, ...
        'MinPeakHeight',   minPeakHeight, ...
        'MinPeakDistance', minPeakDist, ...
        'MinPeakProminence', promRel*(max(f_s)-min(f_s)));
else
    % 简易峰：局部极大点
    loc0 = find(f_s(2:end-1) >= f_s(1:end-2) & f_s(2:end-1) > f_s(3:end)) + 1;
    pk0  = f_s(loc0);
    keep = pk0 >= minPeakHeight;
    loc  = loc0(keep);
    pk   = pk0(keep);
    tmpMin = movmin(f_s, 3, 'Endpoints','shrink');
    prom   = pk - max(min(f_s), tmpMin(loc));   % 简易prominence

end

F.n_peaks = numel(pk);

% 负相段：必须连续至少 Lmin 天才计一次有效翻转
Lmin_neg = 5;
negMask = (f < -eps0);
d = diff([0; negMask; 0]);
neg_st = find(d==1); neg_ed = find(d==-1)-1;
neg_len = neg_ed - neg_st + 1;
F.n_pos_to_neg = sum(neg_len >= Lmin_neg);

% 两峰分波判定：需要相邻两峰之间的谷底足够“深”
isSplitByValley = false;
if numel(loc) >= 2
    for j = 1:(numel(loc)-1)
        seg = f_s(loc(j):loc(j+1));
        vmin = min(seg);
        thresh_valley = max(eps0, min(pk(j), pk(j+1))*(1 - valleyDropRel));
        if vmin <= thresh_valley
            isSplitByValley = true;
            break;
        end
    end
end

% 最终多波标志：严格一些，满足任一严苛条件即可
F.flag_multipart = double( ...
      (F.n_peaks >= 3) ...                       % 峰≥3基本可以认为多波
   || (F.n_peaks >= 2 && isSplitByValley) ...    % 两峰且谷底“够深”
   || (F.n_pos_to_neg >= 2) ...                  % 负相段（持续>=Lmin）出现两次及以上
);

% F.flag_multipart = double( (F.n_peaks >= 2 && isSplitByValley) && (F.n_pos_to_neg >= 2) );

post_len = T - max(F.t_fmax, F.t_inc_slope_max);
F.w_tail_14     = min(1.0, max(0, post_len) / 14);
F.tail_short_14 = double(post_len < 14);
tail_ref_eff    = min(14, round(0.5 * T));
F.w_tail_rel     = min(1.0, max(0, post_len) / max(1, tail_ref_eff));
F.tail_short_rel = double(post_len < tail_ref_eff);

% ---- 生存字段 ----
F.event_time     = iff(isnan(t_cross), T, t_cross);
F.event_observed = double(~isnan(t_cross));
F.event_type     = cross_type;

% ---- meta ----
F.beta = beta; F.gamma = gamma;
end

% ---------- helpers ----------
function v = iff(cond, a, b), if cond, v=a; else, v=b; end, end
function [mx, idx] = safe_max(x)
    if isempty(x), mx = NaN; idx = NaN; return; end
    [mx, idx] = max(x); if isempty(idx), idx = NaN; end
end
function n = count_local_peaks(x)
    x = x(:); N = numel(x); n = 0;
    if N < 3, return; end
    for i = 2:N-1
        if x(i) >= x(i-1) && x(i) > x(i+1), n = n + 1; end
    end
end
