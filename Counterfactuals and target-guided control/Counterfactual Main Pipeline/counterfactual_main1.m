function counterfactual_main1(in_mat, out_root, varargin)
% COUNTERFACTUAL_MAIN  (v3: zero-threshold crossing + timeseries saving)
% 读取 features.mat（需含 Features, f_cell, init_cell）
% 对 f 做反事实编辑（不改步长/长度），统一调用 integrate_seir_odeint_style 求解
% 输出：
%   - <out_root>/cf_results.mat, cf_results.csv   （行：每波×场景）
%   - <out_root>/timeseries/wave_<wave_id>.mat    （含基线与各场景的 Se/S/I/R/inc/f）
%
% 重要口径：以 “0 阈值” 区分正/负相；t_cross 为主峰后首次出现 ≥m_neg 连续 f<0 的起点。

%% ---- 参数 ----
p = inputParser; p.FunctionName = mfilename;
addRequired(p,'in_mat',@(x)ischar(x)||isstring(x));
addRequired(p,'out_root',@(x)ischar(x)||isstring(x));

% 连续点判定（进入负相至少 m_neg 点 <0）
addParameter(p,'m_neg',0, @(x) isnumeric(x) && isscalar(x) && x>=1);
addParameter(p,'m_nz', 0, @(x) isnumeric(x) && isscalar(x) && x>=1);  % 保留不用

% TGC_R2ZERO 控制参数
addParameter(p,'T_h',7,      @(x) isnumeric(x) && isscalar(x) && x>0);     % 目标回落时长
addParameter(p,'auto_T',true, @(x) islogical(x) || ismember(x,[0 1]));      % 短波自适应
addParameter(p,'auto_T_frac',0.6,@(x) isnumeric(x) && isscalar(x) && x>0 && x<=1); % post_len比例
addParameter(p,'t0_mode','t_cross',@(x)ischar(x)||isstring(x));             % 't_cross'|'peak'|数值（未使用）

% 其它
addParameter(p,'save_timeseries',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'verbose',true,      @(x)islogical(x)||ismember(x,[0 1]));
parse(p,in_mat,out_root,varargin{:});
A = p.Results;

%% ---- 读入 ----
if ~isfile(in_mat), error('File not found: %s', in_mat); end
S = load(in_mat,'Features','f_cell','init_cell');
Features  = S.Features;
f_cell    = S.f_cell;
init_cell = S.init_cell;

N = numel(f_cell);
hasGamma = ismember('gamma', Features.Properties.VariableNames);

out_root = char(out_root);
if ~isfolder(out_root), mkdir(out_root); end
ts_dir = fullfile(out_root,'timeseries');
if A.save_timeseries && ~isfolder(ts_dir), mkdir(ts_dir); end

if A.verbose
    fprintf('[counterfactual] Waves=%d → %s\n', N, out_root);
end

%% ---- 场景列表（可按需增删）----
Scenarios = { ...
    struct('name','TGC_R2ZERO',  'type','tgc_r2zero'), ...
    struct('name','TGC_STRONG',  'type','tgc_strong', 'boost',2), ...
    struct('name','NEG_STRONG',  'type','neg_strong', 'scale',2), ...
    struct('name','TGC_LONG',    'type','tgc_long',   'L',3), ...
    struct('name','ADV3',        'type','shift_left', 'k',3), ...
    };

%% ---- 结果表 ----
ResTbl = table();

%% ---- 主循环 ----
for k = 1:N
    try
        % ------- 基础数据 -------
        f0    = ensure_col(double(f_cell{k}));
        init0 = ensure_init(init_cell{k});
        beta  = double(Features.beta(k));
        gamma = hasGamma * double(Features.gamma(k)) + (~hasGamma) * (1/7);
        T     = numel(f0);

        % ------- 0 阈值 t_cross（主峰后，首次 ≥m_neg 连续 f<0）-------
        f_den = movmean(f0, max(1, min(5, 2*floor(0.02*T)+1))); % 轻去抖
        f_den = f0;
        [~, idx_fmax] = max(f_den);
        t_cross = find_t_cross_zero(f_den, idx_fmax, A.m_neg);   % NaN 也可
        t0_pk1 = min(T, idx_fmax);

        % ------- 基线积分（严格用你的积分器） -------
        traj0 = safe_integrate(f0, init0, beta, gamma, 'baseline');
        assert_len(traj0, T, 'baseline');

        tgrid  = traj0.t(:);
        Se0    = traj0.Se(:);
        S0     = traj0.S(:);
        I0     = traj0.I(:);
        Rcomp0 = traj0.R(:);
        inc0   = traj0.I_ode(:);
        C0     = nansum(inc0);

        wave_id  = string(Features.wave_id(k));
        province = string(Features.province(k));
        wave_idx = getfield_if(Features,'wave_idx',k, NaN);

        % ------- timeseries 容器 -------
        if A.save_timeseries
            TS = struct();
            TS.t   = tgrid;
            TS.f   = f0;
            TS.Se  = Se0; TS.S = S0; TS.I = I0; TS.R = Rcomp0;
            TS.inc = inc0;
            TS.meta = struct('beta',beta,'gamma',gamma, ...
                'wave_id',wave_id,'wave_idx',double(wave_idx), ...
                'province',province,'t_cross',double(t_cross),'t0',double(t0_pk1));
        end

        % ------- 各场景 -------
        for s = 1:numel(Scenarios)
            Sc = Scenarios{s};

            % 生成编辑后的 f（0 阈值口径）
            f_cf = apply_scn_f(f0, Sc, t_cross, init0, beta, gamma, A);

            nm = upper(string(Sc.name));
            if startsWith(nm,"TGC")
                % 允许比基线峰值高 10%，且不去尖（避免改变正相形状）
                f_cf = enforce_bounds(f_cf, f0, 'cap_ratio', 1.10, 'despike', false);
            end

            % 再积分（你的积分器；如失败调用内部回退）
            traj  = safe_integrate(f_cf, init0, beta, gamma, sprintf('scenario %s', Sc.name));
            assert_len(traj, T, sprintf('scenario %s', Sc.name));

            inc   = traj.I_ode(:);
            C     = nansum(inc);

            % 结果行
            row = table();
            row.wave_id     = wave_id;
            row.province    = province;
            row.wave_idx    = double(wave_idx);
            row.beta        = beta;
            row.gamma       = gamma;
            row.scenario    = string(Sc.name);
            row.C_base      = C0;
            row.C_cf        = C;
            row.avoided     = C0 - C;
            row.avoided_pct = 100*(C0 - C)/max(C0, eps);
            row.t_cross     = double(t_cross);

            ResTbl = [ResTbl; row]; %#ok<AGROW>

            % 保存时序
            if A.save_timeseries
                tag = char(Sc.name);
                TS.("f_"+tag)   = f_cf;
                TS.("Se_"+tag)  = traj.Se(:);
                TS.("S_"+tag)   = traj.S(:);
                TS.("I_"+tag)   = traj.I(:);
                TS.("R_"+tag)   = traj.R(:);
                TS.("inc_"+tag) = inc;
            end
        end

        if A.save_timeseries
            save(fullfile(ts_dir, sprintf('wave_%s.mat', sanitize(char(wave_id)))), 'TS');
        end

    catch ME
        warning('Wave %d (%s) failed: %s @ %s:%d', k, string(Features.wave_id(k)), ...
            ME.message, topfile(ME), topline(ME));
    end
end

%% ---- 保存索引 ----
Results = ResTbl;
save(fullfile(out_root,'cf_results.mat'), 'Results','Scenarios');
writetable(ResTbl, fullfile(out_root,'cf_results.csv'));

if A.verbose
    fprintf('[counterfactual] Done. Waves=%d, scenarios=%d → %s\n', N, numel(Scenarios), out_root);
end
end

%% ================= 工具函数 =================
function v = ensure_col(v), v = v(:); end
function s = sanitize(s), s = regexprep(s, '[^a-zA-Z0-9-_]', '_'); end
function f = topfile(ME), if isempty(ME.stack), f='?'; else, f=ME.stack(1).name; end, end
function l = topline(ME), if isempty(ME.stack), l=0; else, l=ME.stack(1).line; end, end

function val = getfield_if(T, name, idx, fallback)
if ismember(name, T.Properties.VariableNames), val = T.(name)(idx);
else, val = fallback; end
end

function init = ensure_init(x)
% 支持 struct(Se0,S0,I0,R0)/别名 或 4x1 向量
if isstruct(x)
    init = struct('Se0',0,'S0',0,'I0',0,'R0',0);
    if isfield(x,'Se0'), init.Se0 = double(x.Se0);
    elseif isfield(x,'E0'), init.Se0 = double(x.E0);
    elseif isfield(x,'Se'), init.Se0 = double(x.Se);
    end
    if isfield(x,'S0'), init.S0 = double(x.S0); elseif isfield(x,'S'), init.S0 = double(x.S); end
    if isfield(x,'I0'), init.I0 = double(x.I0); elseif isfield(x,'I'), init.I0 = double(x.I); end
    if isfield(x,'R0'), init.R0 = double(x.R0); elseif isfield(x,'R'), init.R0 = double(x.R); end
else
    v = double(x(:)); v(end+1:4,1) = 0;
    init = struct('Se0',v(1),'S0',v(2),'I0',v(3),'R0',v(4));
end
end

function assert_len(traj, T, tag)
ok = (numel(traj.t)==T) && (numel(traj.S)==T) && (numel(traj.I)==T) && (numel(traj.R)==T) && (numel(traj.I_ode)==T);
if ~ok
    error('Length mismatch (%s): t=%d S=%d I=%d R=%d I_ode=%d, expected %d', ...
        tag, numel(traj.t), numel(traj.S), numel(traj.I), numel(traj.R), numel(traj.I_ode), T);
end
end

function t_cross = find_t_cross_zero(f, idx_fmax, m_neg)
% 主峰之后，第一次出现 >=m_neg 个连续 f<0 的起点（阈值=0）
T = numel(f); t_cross = NaN;
if T > idx_fmax
    run = 0;
    for t = (idx_fmax+1):T
        if f(t) < 0
            run = run + 1;
            if run >= m_neg
                t_cross = t - m_neg;
                return;
            end
        else
            run = 0;
        end
    end
end
end

function y = getfield_default(S, name, val)
if isfield(S, name), y = S.(name); else, y = val; end
end

%% === 场景编辑（基于 0 阈值） ===
function f_cf = apply_scn_f(f0, Sc, t_cross, init0, beta, gamma, A)
% 单输出，保持长度/步长不变

T  = numel(f0); f0 = f0(:); f_cf = f0;
% —— 与主循环一致的轻去抖（窗长 ≈ max(1, min(5, 2*floor(0.02*T)+1))）
win = max(1, min(5, 2*floor(0.02*T)+1));
f_s = movmean(f0, win);


switch upper(string(Sc.name))
    case "BASELINE"
        % 不变

    case "ADV3"
        k = getfield_default(Sc,'k',3);
        k = min(T-1, max(0, round(k)));
        f_cf = [f0(1+k:end); zeros(k,1)];

    case "NEG_STRONG"
        scale = getfield_default(Sc,'scale',2);
        mneg  = (f0 < 0);
        f_cf  = f0;
        f_cf(mneg) = scale * f0(mneg);

    case "TGC_R2ZERO"
        % 起控点：峰后首次 f<=0；若未到负相，则用峰值
        t0   = pick_t0_zero(f_s);
        Th   = pick_Th(A, t0, T);             % 你已有的短波自适应函数
        prm  = struct('T_h',Th,'t0',t0,'lambda',0.55, ...
            'max_iter',3,'f_bounds',[-3,3], ...
            'deepen_neg',0.00);      % 不触碰负相
        f_cf = synthesize_tgc_r2zero_masked_zero(f0, init0, beta, gamma, prm);

    case "TGC_STRONG"
        t0   = pick_t0_zero(f_s);
        Th   = max(5, round(0.7 * pick_Th(A, t0, T))); % 更快回零
        prm  = struct('T_h',Th,'t0',t0,'lambda',0.85, ...
            'max_iter',4,'f_bounds',[-3,3], ...
            'deepen_neg',0.10);      % 仅“更负一点”（10%），绝不往 0 拉
        f_cf = synthesize_tgc_r2zero_masked_zero(f0, init0, beta, gamma, prm);

    case "TGC_LONG"
        % 让负相更持久：稍早起控，但仍只压正相
        t0   = max(1, pick_t0_zero(f_s) - 2);
        Th   = max(7, round(0.8 * pick_Th(A, t0, T)));
        prm  = struct('T_h',Th,'t0',t0,'lambda',0.75, ...
            'max_iter',3,'f_bounds',[-3,3], ...
            'deepen_neg',0.00);
        f_cf = synthesize_tgc_r2zero_masked_zero(f0, init0, beta, gamma, prm);


    otherwise
        % 未知场景：保持不变
end

end

function t0 = pick_t0_zero(f0)
% 峰后首次 f<=0 的时刻；若全程非负，就取峰值处
T = numel(f0);
[~, tpk] = max(f0);
t0 = tpk;
for tt = tpk:T
    if f0(tt) <= 0
        t0 = tt; break;
    end
end

t0 = tpk+1;

end



function Th = pick_Th(A, t0, T)
post_len = max(0, T - t0 + 1);
if isfield(A,'auto_T') && A.auto_T
    Th = min(A.T_h, max(3, round(A.auto_T_frac * max(1,post_len))));
else
    Th = A.T_h;
end
end


%% === 公式驱动 TGC：让 r(t) 在 T_h 内回到 0（逐次替代法） ===
function f_cf = synthesize_tgc_r2zero_masked_zero(f_base, init0, beta, gamma, prm)
% 目标导向控制（零阈值版）：仅在 t>=t0 的“正相”(f>0) 位置更新 f，使 r(t) 在 T_h 内
% 指数回落到 ~0；对 f<=0 的位置不动（Strong 可选择加深负相一点点）。

% ------- 参数 -------
T       = numel(f_base);  f_base = f_base(:);
t0      = max(1, min(T, round(prm.t0)));
Th      = max(3, round(prm.T_h));
lambda  = min(max(prm.lambda,0),1);
lohi    = prm.f_bounds; lo = lohi(1); hi = lohi(2);
maxIter = max(1, round(prm.max_iter));
deepNeg = max(0, getfield_default(prm,'deepen_neg',0.0));   % e.g., 0.1 => 负相×1.1
k       = 1 / Th;

% ------- 初始 ------
f_curr = f_base;

for it = 1:maxIter
    % 1) 用当前 f 积分，得到 r,s,i
    traj = integrate_seir_odeint_style(f_curr, init0, beta, gamma);
    if numel(traj.t) ~= T
        error('tgc_r2zero_zero:early_stop','integrator stopped at t=%d (expect %d)', numel(traj.t), T);
    end
    S = traj.S(:); I = traj.I(:); R = traj.R(:);
    N = S + I + R; N(N<=0)=eps;
    s = S./N; i = I./N;
    r = beta.*s - gamma;

    % 2) 目标 r_tar（只定义下降轨，但等下只在“正相掩码”处使用）
    r_tar = r;
    if t0 <= T
        idx = (t0:T)';
        r_tar(idx) = r(t0) .* exp(-k*(idx - t0));
    end

    % 3) 反推 f*（全长计算，但只对掩码点更新）
    denom  = max(1e-8, 1 - s);
    f_star = (N ./ denom) .* ( - (k/beta).*r_tar + beta .* s .* i );

    % 4) 掩码与更新：只更新“正相” -> f_base>0 且 t>=t0
    mask_update = false(T,1);
    if t0 <= T
        mask_update(t0:T) = (f_base(t0:T) > 0);
    end

    f_new = f_curr;

    if any(mask_update)
        seg = f_curr(mask_update);
        seg_target = f_star(mask_update);
        seg = lambda.*seg_target + (1-lambda).*seg;
        seg = min(max(seg, lo), hi);
        % 保护：不允许把“本该正相”的点改成负值（更符合“只压不翻负”的直觉）
        seg = max(seg, 0);
        f_new(mask_update) = seg;
    end

    % 5) 负相点：完全不动；Strong 可“加深一点”（更负）
    if deepNeg > 0 && t0 <= T
        mneg = (f_base<=0);
        mneg(1:t0-1) = false;      % 仅 t>=t0
        if any(mneg)
            f_new(mneg) = min(f_new(mneg), (1+deepNeg)*f_base(mneg));
        end
    end

    % 6) 不做平滑
    f_curr = f_new;
end

f_cf = f_curr;
end



%% === 积分（严格→回退） ===
function traj = safe_integrate(f, init0, beta, gamma, tag)
T = numel(f);
try
    traj = integrate_seir_odeint_style(f, init0, beta, gamma);
    if numel(traj.t) ~= T
        error('safe_integrate:early_stop','early stop at t=%d (expect %d)', numel(traj.t), T);
    end
catch
    % 回退：放宽容差（保持整日输出）
    tspan = 0:(T-1);
    y0    = [init0.Se0; init0.S0; init0.I0; init0.R0];

    opts2 = odeset('RelTol',1e-7,'AbsTol',1e-10);  % 略放宽
    [tt, YY] = ode45(@rhs, tspan, y0, opts2);

    if numel(tt) ~= T
        error('safe_integrate:still_early_stop','t=%d, expect=%d', numel(tt), T);
    end

    Se = YY(:,1); S = YY(:,2); I = YY(:,3); R = YY(:,4);
    den = S + I + R; den(den<=0) = eps;
    I_ode = beta .* S .* I ./ den;

    traj = struct('t',tt,'Se',Se,'S',S,'I',I,'R',R,'I_ode',I_ode);
    warning('[safe_integrate] relaxed options used for %s', tag);
end

    function dydt = rhs(tt, y)
        Se = y(1); S = y(2); I = y(3); R = y(4);
        denom = S + I + R; if denom <= 0, denom = eps; end
        inc = beta * S * I / denom;
        ft  = NN_f(tt);
        dydt = [-ft;        % dSe/dt
            ft - inc;   % dS/dt
            inc - gamma*I;  % dI/dt
            gamma*I];       % dR/dt
    end

    function val = NN_f(tt)
        rr  = round_half_to_even(tt);
        rr  = max(0, min(T-1, rr));
        idx = rr + 1;
        val = f(idx);
    end
end

function r = round_half_to_even(x)
fx   = floor(x);
frac = x - fx;
r    = fx;
gt   = frac > 0.5;
eq   = abs(frac - 0.5) <= 1e-12;
r(gt) = fx(gt) + 1;
odd   = mod(fx,2) ~= 0;
bump  = eq & odd;
r(bump) = fx(bump) + 1;
end

function f2 = enforce_bounds(f, f_ref, varargin)
% ENFORCE_BOUNDS  可选的温和兜底：限幅/去尖都默认关闭
p = inputParser; p.FunctionName = mfilename;
addParameter(p,'cap_ratio', inf, @(x)isnumeric(x)&&isscalar(x)&&x>0);   % =Inf 不限幅
addParameter(p,'despike',   false, @(x)islogical(x)||ismember(x,[0 1]));
parse(p,varargin{:});
A = p.Results;

f = f(:); f_ref = f_ref(:); T = numel(f);

% —— 仅当 cap_ratio < Inf 时才限幅（对称）——
if isfinite(A.cap_ratio)
    mu   = mean(f_ref);
    madm = mean(abs(f_ref - mu));
    fmax = max(abs(f_ref));
    cap0 = max(3*madm,  A.cap_ratio * fmax);   % 比如 cap_ratio=1.1 会允许比基线峰高 10%
    f    = min(max(f, -cap0), cap0);
end

% —— 仅当 despike=true 时才去尖 —— 
if A.despike && T>=3
    jump = [0; diff(f)];
    [~, i1] = max(abs(jump));
    i1 = max(2, min(T-1, i1));
    f(i1) = 0.5*(f(i1-1) + f(i1+1));
end

f2 = f;
end

