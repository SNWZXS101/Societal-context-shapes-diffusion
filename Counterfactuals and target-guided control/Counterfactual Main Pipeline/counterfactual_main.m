function counterfactual_main(in_mat, out_root, varargin)
% COUNTERFACTUAL_MAIN  (v4: zero-threshold, adaptive mask, neg floor, robust integrator)
% 读取 features.mat（需含 Features, f_cell, init_cell）
% 对 f 做反事实编辑（不改步长/长度），统一调用 integrate_seir_odeint_style 求解
% 输出：
%   - <out_root>/cf_results.mat, cf_results.csv   （行：每波×场景）
%   - <out_root>/timeseries/wave_<wave_id>.mat    （含基线与各场景的 Se/S/I/R/inc/f）
%
% 重要口径：
%   • 正/负相阈值固定为 0
%   • t0 选取：主峰后第 1 天（在去抖后的 f 上找主峰，再 +1；超界则夹到 T）
%   • TGC 允许 f 适度为负，但不低于“负向地板”（默认 -f_peak，可配）
%   • 积分：严格 → ode45(宽容) → ode15s(刚性) 三段回退；始终输出整日网格

%% ---- 参数 ----
p = inputParser; p.FunctionName = mfilename;
addRequired(p,'in_mat',@(x)ischar(x)||isstring(x));
addRequired(p,'out_root',@(x)ischar(x)||isstring(x));

% 连续点判定（进入负相至少 m_neg 点 <0）——仅用于元数据 t_cross，可不依赖
addParameter(p,'m_neg',1, @(x) isnumeric(x) && isscalar(x) && x>=1);

% TGC 控制参数
addParameter(p,'T_h',10,      @(x) isnumeric(x) && isscalar(x) && x>0);      % 目标回落时长
addParameter(p,'auto_T',true, @(x) islogical(x) || ismember(x,[0 1]));       % 短波自适应
addParameter(p,'auto_T_frac',0.6,@(x) isnumeric(x) && isscalar(x) && x>0 && x<=1);

% TGC：允许负向的“地板”与更新掩码
addParameter(p,'neg_floor_ratio',1.0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);  % 地板 = -ratio * f_peak
addParameter(p,'neg_floor_abs',[],    @(x)isempty(x) || (isnumeric(x)&&isscalar(x))); % 绝对地板（优先）
addParameter(p,'update_mask_mode','adaptive_floor', @(x)ischar(x)||isstring(x));
% 选项：'adaptive_floor'(默认) | 'all_t_ge_t0' | 'baseline_pos'(旧口径)

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
    struct('name','TGC_STRONG',  'type','tgc_strong'), ...
    struct('name','TGC_LONG',    'type','tgc_long',   'L',3), ...
    struct('name','NEG_STRONG',  'type','neg_strong', 'scale',1.5), ...
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

        % ------- 去抖 & t_cross（仅元数据用） -------
        win    = max(1, min(5, 2*floor(0.02*T)+1));   % 约 2%T 的奇窗
        f_s    = movmean(f0, win);
        [~, idx_fmax_den] = max(f_s);
        t_cross = find_t_cross_zero(f_s, idx_fmax_den, A.m_neg);   % 可能为 NaN

        % ------- t0：主峰后 +1（在去抖后的 f 上找峰） -------
        t0 = pick_t0_pk1(f_s);   % = min(T, tpk+1)

        % ------- 基线积分（严格 → 宽容回退） -------
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
                'province',province,'t_cross',double(t_cross), ...
                't0_pk1',double(t0));
        end

        % ------- 各场景 -------
        for s = 1:numel(Scenarios)
            Sc = Scenarios{s};

            % 生成编辑后的 f（0 阈值口径 + t0=峰后+1）
            f_cf = apply_scn_f(f0, f_s, Sc, t0, init0, beta, gamma, A);

            % TGC：允许比基线峰值高 10%，且不去尖（避免改变正相形状）
            nm = upper(string(Sc.name));
            if startsWith(nm,"TGC")
                f_cf = enforce_bounds(f_cf, f0, 'cap_ratio', 1.10, 'despike', false);
            end

            % 再积分
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

                % 写入该场景实际使用的 t0/Th（如适用）
                if startsWith(nm,"TGC")
                    % 将 pick_Th 的结果也记入 meta
                    Th_used = compute_Th_used(Sc, A, t0, T);
                    TS.meta.("t0_"+tag) = double(t0);
                    TS.meta.("Th_"+tag) = double(Th_used);
                end
            end
        end

        if A.save_timeseries
            save(fullfile(ts_dir, sprintf('wave_%s.mat', sanitize(char(wave_id)))), 'TS');
        end

    catch ME
        warning('[counterfactual] Wave %d (%s) failed: %s @ %s:%d', ...
            k, string(Features.wave_id(k)), ME.message, topfile(ME), topline(ME));
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
if istable(T)
    if ismember(name, T.Properties.VariableNames), val = T.(name)(idx);
    else, val = fallback; end
else
    val = fallback;
end
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
                t_cross = t - m_neg + 1;   % 起点
                return;
            end
        else
            run = 0;
        end
    end
end
end

function t0 = pick_t0_pk1(f_s)
% 峰后 +1（用去抖 f_s 找峰）
T = numel(f_s);
[~, tpk] = max(f_s);
t0 = min(T, tpk + 1);
end

function Th = pick_Th(A, t0, T)
% auto_T=false 时严格使用传入的 T_h；开启自适应时再截断
post_len = max(0, T - t0 + 1);
if isfield(A,'auto_T') && ~A.auto_T
    Th = A.T_h;
else
    Th = min(A.T_h, max(3, round(A.auto_T_frac * max(1,post_len))));
end
end

function Th_used = compute_Th_used(Sc, A, t0, T)
% 用于写 meta：给 TGC 家族返回实际用的 Th
nm = upper(string(Sc.name));
if startsWith(nm,"TGC")
    Th_used = pick_Th(A, t0, T);
    if nm=="TGC_STRONG"
        Th_used = max(5, round(0.7*Th_used));
    elseif nm=="TGC_LONG"
        Th_used = max(7, round(0.8*Th_used));
    end
else
    Th_used = NaN;
end
end

%% === 场景编辑（0 阈值口径 + t0=峰后+1） ===
function f_cf = apply_scn_f(f0, f_s, Sc, t0, init0, beta, gamma, A)
% f0: 原始；f_s: 去抖后（仅用于 t0 推断；已在主循环完成）
T  = numel(f0);
f0 = f0(:);
f_cf = f0;    % 默认不变

switch upper(string(Sc.name))
    case "BASELINE"
        % 不变

    case "ADV3"
        k = getfield_default(Sc,'k',3);
        k = min(T-1, max(0, round(k)));
        f_cf = [f0(1+k:end); zeros(k,1)];

    case "NEG_STRONG"
        scale = getfield_default(Sc,'scale',1.5);
        mneg  = (f0 < 0);
        f_cf  = f0;
        f_cf(mneg) = scale * f0(mneg);

    case "TGC_R2ZERO"
        Th   = pick_Th(A, t0, T);
        prm  = struct( ...
            'T_h',Th,'t0',t0,'lambda',0.55, ...
            'max_iter',3,'f_bounds',[-3,3], ...
            'deepen_neg',0.00, ...
            'neg_floor_ratio', A.neg_floor_ratio, ...
            'neg_floor_abs',   A.neg_floor_abs, ...
            'update_mask_mode',string(A.update_mask_mode));
        f_cf = synthesize_tgc_r2zero_masked_zero(f0, init0, beta, gamma, prm);

    case "TGC_STRONG"
        Th0  = pick_Th(A, t0, T);
        Th   = max(5, round(0.7 * Th0));     % 更快回零
        prm  = struct( ...
            'T_h',Th,'t0',t0,'lambda',0.85, ...
            'max_iter',4,'f_bounds',[-3,3], ...
            'deepen_neg',0.10, ...
            'neg_floor_ratio', A.neg_floor_ratio, ...
            'neg_floor_abs',   A.neg_floor_abs, ...
            'update_mask_mode',string(A.update_mask_mode));
        f_cf = synthesize_tgc_r2zero_masked_zero(f0, init0, beta, gamma, prm);

    case "TGC_LONG"
        Th0  = pick_Th(A, t0, T);
        t0L  = max(1, t0 - 2);
        Th   = max(7, round(0.8 * Th0));
        prm  = struct( ...
            'T_h',Th,'t0',t0L,'lambda',0.75, ...
            'max_iter',3,'f_bounds',[-3,3], ...
            'deepen_neg',0.00, ...
            'neg_floor_ratio', A.neg_floor_ratio, ...
            'neg_floor_abs',   A.neg_floor_abs, ...
            'update_mask_mode',string(A.update_mask_mode));
        f_cf = synthesize_tgc_r2zero_masked_zero(f0, init0, beta, gamma, prm);

    otherwise
        % 未知场景：保持不变
end

end

%% === 公式驱动 TGC（零阈值版，允许负到地板） ===
function f_cf = synthesize_tgc_r2zero_masked_zero(f_base, init0, beta, gamma, prm)
% t>=t0 区间内朝目标 r(t) 指数回零更新 f；
% —— 自适应更新掩码：依据每次迭代的 f_curr（非基线）；
% —— 允许为负，但不低于“负向地板”（相对峰值或绝对值）；
% —— 负相“更负一点”仍由 deepen_neg 控制。

% ------- 参数 -------
T       = numel(f_base);  f_base = f_base(:);
t0      = max(1, min(T, round(prm.t0)));
Th      = max(3, round(prm.T_h));
lambda  = min(max(getfield_default(prm,'lambda',0.55),0),1);
lohi    = getfield_default(prm,'f_bounds',[-3,3]); lo = lohi(1); hi = lohi(2);
maxIter = max(1, round(getfield_default(prm,'max_iter',3)));
deepNeg = max(0,   getfield_default(prm,'deepen_neg',0.0));
maskMode= lower(string(getfield_default(prm,'update_mask_mode','adaptive_floor')));

% ------- 负向地板（绝对优先于比例） -------
if isfield(prm,'neg_floor_abs') && ~isempty(prm.neg_floor_abs)
    neg_floor = double(prm.neg_floor_abs);
else
    ratio     = max(0, getfield_default(prm,'neg_floor_ratio',1.0));
    fpos_pk   = max(0, max(f_base));         % 正相峰值
    neg_floor = -ratio * fpos_pk;            % 例如 -1*f_peak
end
neg_floor = max(neg_floor, lo);              % 不越过全局下界
k      = 1/Th;
f_curr = f_base;

for it = 1:maxIter
    % ---- 积分，得到 r,s,i ----
    traj = integrate_seir_odeint_style(f_curr, init0, beta, gamma);
    if numel(traj.t) ~= T
        error('tgc_r2zero_zero:early_stop','integrator stopped at t=%d (expect %d)', numel(traj.t), T);
    end
    S = traj.S(:); I = traj.I(:); R = traj.R(:);
    N = S + I + R; N(N<=0)=eps;
    s = S./N; i = I./N;
    r = beta.*s - gamma;

    % ---- 目标 r_tar：t>=t0 指数回零 ----
    r_tar = r;
    if t0 <= T
        idx = (t0:T)'; r_tar(idx) = r(t0) .* exp(-k*(idx - t0));
    end

    % ---- 反推 f* ----
    denom  = max(1e-8, 1 - s);
    f_star = (N ./ denom) .* ( - (k/beta).*r_tar + beta .* s .* i );

    % ---- 更新掩码（自适应/全开/旧口径）----
    mask_update = false(T,1);
    if t0 <= T
        switch maskMode
            case 'baseline_pos'
                mask_update(t0:T) = (f_base(t0:T) > 0);
            case 'all_t_ge_t0'
                mask_update(t0:T) = true;
            otherwise % 'adaptive_floor'
                mask_update(t0:T) = (f_curr(t0:T) >= (neg_floor - 1e-9));
        end
    end

    % ---- 朝目标推进 + 全局限幅 + 负向地板（无 >=0 硬夹）----
    f_new = f_curr;
    if any(mask_update)
        seg        = f_curr(mask_update);
        seg_target = f_star(mask_update);

        seg = lambda.*seg_target + (1-lambda).*seg;
        seg = min(max(seg, lo), hi);   % 全局边界（如 [-3,3]）
        seg = max(seg, neg_floor);     % ★ 允许负，但不低于地板

        f_new(mask_update) = seg;
    end

    % ---- 负相可“更负一点”（Strong 场景）----
    if deepNeg > 0 && t0 <= T
        mneg = (f_base<=0); mneg(1:t0-1) = false;
        if any(mneg)
            f_new(mneg) = min(f_new(mneg), (1+deepNeg)*f_base(mneg));
        end
    end

    f_curr = f_new;
end

f_cf = f_curr;
end

%% === 积分（严格→ode45宽容→ode15s刚性） ===
function traj = safe_integrate(f, init0, beta, gamma, tag)
T = numel(f);

% 1) 严格：项目口径的自定义积分器
try
    traj = integrate_seir_odeint_style(f, init0, beta, gamma);
    if numel(traj.t) ~= T
        error('safe_integrate:early_stop','early stop at t=%d (expect %d)', numel(traj.t), T);
    end
    return;
catch
    % fall through
end

% 2) 宽容：ode45 + 非负裁剪右端（整日输出）
try
    tspan = 0:(T-1);
    y0    = [init0.Se0; init0.S0; init0.I0; init0.R0];
    opts2 = odeset('RelTol',5e-6,'AbsTol',5e-9);  % 略放宽
    [tt, YY] = ode45(@rhs_relaxed, tspan, y0, opts2);
    if numel(tt) ~= T, error('len'); end
    Se=YY(:,1); S=YY(:,2); I=YY(:,3); R=YY(:,4);
    den = S+I+R; den(den<=0)=eps; I_ode = beta.*S.*I./den;
    traj = struct('t',tt,'Se',Se,'S',S,'I',I,'R',R,'I_ode',I_ode);
    warning('[safe_integrate] relaxed ode45 used for %s', tag);
    return;
catch
    % fall through
end

% 3) 刚性：ode15s + NonNegative（整日输出）
tspan = 0:(T-1);
y0    = [init0.Se0; init0.S0; init0.I0; init0.R0];
opts3 = odeset('RelTol',1e-7,'AbsTol',1e-10,'NonNegative',1:4);
[tt, YY] = ode15s(@rhs_relaxed, tspan, y0, opts3);
if numel(tt) ~= T
    error('safe_integrate:still_early_stop','t=%d, expect=%d', numel(tt), T);
end
Se=YY(:,1); S=YY(:,2); I=YY(:,3); R=YY(:,4);
den = S+I+R; den(den<=0)=eps; I_ode = beta.*S.*I./den;
traj = struct('t',tt,'Se',Se,'S',S,'I',I,'R',R,'I_ode',I_ode);
warning('[safe_integrate] ode15s fallback used for %s', tag);

    function dydt = rhs_relaxed(tt, y)
        % 与 integrate_seir_odeint_style 的右端一致，但对状态做非负裁剪防抖
        Se = max(0.0, y(1)); S = max(0.0, y(2));
        I  = max(0.0, y(3)); R = max(0.0, y(4));
        den = S + I + R; if den <= 0, den = eps; end
        inc = beta * S * I / den;
        ft  = NN_f(tt);
        dydt = [-ft; ft - inc; inc - gamma*I; gamma*I];
    end

    function val = NN_f(tt)
        rr  = round_half_to_even(tt);
        rr  = max(0, min(T-1, rr));
        idx = rr + 1;
        val = f(idx);
    end
end

%% === 项目口径的“严格积分器” ===
function out = integrate_seir_odeint_style(f, init, beta, gamma)
% t = 0:1:T-1；f 采用最近邻（ties-to-even）取样；状态参与传播项前做非负裁剪
T = numel(f);
tspan = 0:(T-1);
y0 = [max(0,init.Se0); max(0,init.S0); max(0,init.I0); max(0,init.R0)];
useNonneg = true;
    function idx = idx_from_t(tt)
        rr  = round_half_to_even(tt);
        rr  = max(0, min(T-1, rr));
        idx = rr + 1;
    end
    function val = NN_f(tt), val = f(idx_from_t(tt)); end
    function dydt = rhs(tt, y)
        Se = max(0.0, y(1)); S = max(0.0, y(2));
        I  = max(0.0, y(3)); R = max(0.0, y(4));
        den = S + I + R; if den <= 0, den = eps; end
        inc = beta * S * I / den;
        ft  = NN_f(tt);
        dydt = [-ft; ft - inc; inc - gamma*I; gamma*I];
    end
opts = odeset('RelTol',1e-7,'AbsTol',1e-10);
if useNonneg, opts = odeset(opts,'NonNegative',1:4); end
[tt, YY] = ode45(@rhs, tspan, y0, opts);
Se = max(0, YY(:,1)); S = max(0, YY(:,2));
I  = max(0, YY(:,3)); R = max(0, YY(:,4));
den = S + I + R; den(den<=0) = eps;
I_ode = beta .* S .* I ./ den;
out = struct('t',tt, 'Se',Se, 'S',S, 'I',I, 'R',R, 'I_ode',I_ode);
end

%% === 辅助 ===
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
% 可选兜底：限幅/去尖默认关闭（TGC 可设 cap_ratio=1.10, despike=false）
p = inputParser; p.FunctionName = mfilename;
addParameter(p,'cap_ratio', inf, @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'despike',   false, @(x)islogical(x)||ismember(x,[0 1]));
parse(p,varargin{:});
A = p.Results;
f = f(:); f_ref = f_ref(:); T = numel(f);

if isfinite(A.cap_ratio)
    mu   = mean(f_ref);
    madm = mean(abs(f_ref - mu));
    fmax = max(abs(f_ref));
    cap0 = max(3*madm,  A.cap_ratio * fmax);
    f    = min(max(f, -cap0), cap0);
end

if A.despike && T>=3
    jump = [0; diff(f)];
    [~, i1] = max(abs(jump));
    i1 = max(2, min(T-1, i1));
    f(i1) = 0.5*(f(i1-1) + f(i1+1));
end
f2 = f;
end

function y = getfield_default(S, name, val)
if isstruct(S) && isfield(S, name), y = S.(name); else, y = val; end
end
