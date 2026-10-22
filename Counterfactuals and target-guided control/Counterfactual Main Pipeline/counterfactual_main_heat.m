function counterfactual_main_heat(in_mat, out_root, varargin)
% COUNTERFACTUAL_MAIN_HEAT
% 扫 NEG(eta) 与 TGC(t0+Δ) 的 %Avoid 热图所需数据（不保存时序）
% 口径：
%   - 基线 t0 := 主峰后 +1 天（t_pk+1），TGC 的 Th 固定（默认 10；可 auto_T）
%   - NEG: 仅放大/缩小负相（f<0）部分：f_cf = (1+eta)*f_neg + f_pos，eta∈[-1,2]
%
% 输出：
%   <out_root>/heat_results.mat  （结构 Heat，含 .NEG 与 .TGC_T0 两张长表）
%   <out_root>/heat_NEG.csv, <out_root>/heat_TGC_T0.csv

% ---------- 参数 ----------
p = inputParser; p.FunctionName = mfilename;
addRequired(p,'in_mat',@(x)ischar(x)||isstring(x));
addRequired(p,'out_root',@(x)ischar(x)||isstring(x));
addParameter(p,'EtaMin', -1,  @isnumeric);
addParameter(p,'EtaMax',  2,  @isnumeric);
addParameter(p,'EtaStep', 0.5,@isnumeric);
addParameter(p,'T0DeltaMin', -3, @isnumeric);
addParameter(p,'T0DeltaMax',  3, @isnumeric);
addParameter(p,'T0DeltaStep', 1, @isnumeric);

% TGC 控制固定口径（可按需改）
addParameter(p,'T_h',10,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'auto_T',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'auto_T_frac',0.6,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<=1);

% 允许负相（相对峰值的地板；例如 0.5 表示不低于 -0.5*|f_peak|）
addParameter(p,'neg_floor_rel',0.50,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'deepNeg',0.10,@(x)isnumeric(x)&&isscalar(x)&&x>=0);  % 仅在负相可选“稍加深”

addParameter(p,'Verbose',true,@(x)islogical(x)||ismember(x,[0 1]));
parse(p,in_mat,out_root,varargin{:});
A = p.Results;

in_mat  = char(in_mat);
out_root= char(out_root);
if ~isfile(in_mat), error('File not found: %s', in_mat); end
if ~isfolder(out_root), mkdir(out_root); end

% ---------- 读取 ----------
S = load(in_mat,'Features','f_cell','init_cell');
Features  = S.Features;
f_cell    = S.f_cell;
init_cell = S.init_cell;
N = numel(f_cell);
hasGamma = ismember('gamma', Features.Properties.VariableNames);

% 解析分组标签（variant / cluster）
[variant, cluster] = resolve_labels(Features);

% 参数网格
eta_list   = A.EtaMin:A.EtaStep:A.EtaMax;
dt0_list   = A.T0DeltaMin:A.T0DeltaStep:A.T0DeltaMax;

% 预分配结果表（长表）
T_NEG = table('Size',[0 7], 'VariableTypes',{'string','string','double','double','double','double','double'}, ...
    'VariableNames',{'wave_id','province','cluster','eta','avoided_pct','beta','gamma'});
T_T0  = table('Size',[0 8], 'VariableTypes',{'string','string','double','double','double','double','double','double'}, ...
    'VariableNames',{'wave_id','province','cluster','dt0','avoided_pct','beta','gamma','Th_used'});

if A.Verbose
    fprintf('[cf_heat] Waves=%d, NEG |eta|=%d, TGC_T0 |d|=%d → %s\n', ...
        N, numel(eta_list), numel(dt0_list), out_root);
end

% ---------- 主循环 ----------
for k = 1:N
    try
        f0    = ensure_col(double(f_cell{k}));
        init0 = ensure_init(init_cell{k});
        beta  = double(Features.beta(k));
        gamma = hasGamma * double(Features.gamma(k)) + (~hasGamma) * (1/7);
        T     = numel(f0);

        % 主峰+1 的基线 t0
        [~,tpk] = max(f0);
        t0_pk1  = min(T, tpk+1);

        % 基线积分
        traj0 = safe_integrate(f0, init0, beta, gamma, 'baseline');
        assert_len(traj0, T, 'baseline');
        C0 = nansum(traj0.I_ode);

        wave_id  = string(Features.wave_id(k));
        province = string(Features.province(k));
        clu_k    = cluster(k);

        % -------- NEG 扫 eta --------
        for e = 1:numel(eta_list)
            eta = eta_list(e);
            f_cf = apply_neg_eta(f0, eta);
            traj = safe_integrate(f_cf, init0, beta, gamma, 'NEG');
            assert_len(traj, T, 'NEG');
            C = nansum(traj.I_ode);
            row = {wave_id, province, double(clu_k), double(eta), 100*(C0-C)/max(C0,eps), beta, gamma};
            T_NEG = [T_NEG; row]; %#ok<AGROW>
        end

        % -------- TGC 扫 Δt0 --------
        for d = 1:numel(dt0_list)
            dt0   = dt0_list(d);
            t0    = max(1, min(T, t0_pk1 + round(dt0)));
            Th    = pick_Th_struct(struct('T_h',A.T_h,'auto_T',A.auto_T,'auto_T_frac',A.auto_T_frac), t0, T);
            prm   = struct('T_h',Th,'t0',t0,'lambda',0.70, ...
                           'f_bounds',[-Inf,Inf], ...  % 不做数值限幅；稳定性交给积分器+地板
                           'neg_floor_rel',A.neg_floor_rel,'deepNeg',A.deepNeg);
            f_cf  = synthesize_tgc_r2zero_masked_zero(f0, init0, beta, gamma, prm);
            traj  = safe_integrate(f_cf, init0, beta, gamma, 'TGC_T0');
            assert_len(traj, T, 'TGC_T0');
            C     = nansum(traj.I_ode);
            row   = {wave_id, province, double(clu_k), double(dt0), 100*(C0-C)/max(C0,eps), beta, gamma, double(Th)};
            T_T0  = [T_T0; row]; %#ok<AGROW>
        end

    catch ME
        warning('[cf_heat] Wave %d (%s) failed: %s @ %s:%d', k, string(Features.wave_id(k)), ...
            ME.message, topfile(ME), topline(ME));
    end
end

% ---------- 保存 ----------
Heat = struct();
Heat.NEG = T_NEG;
Heat.TGC_T0 = T_T0;
Heat.Param = struct('eta_list',eta_list,'dt0_list',dt0_list);
Heat.Meta  = struct('variant',{variant},'cluster',cluster);

save(fullfile(out_root,'heat_results.mat'), 'Heat','-v7.3');
writetable(T_NEG, fullfile(out_root,'heat_NEG.csv'));
writetable(T_T0,  fullfile(out_root,'heat_TGC_T0.csv'));

if A.Verbose
    fprintf('[cf_heat] Saved → %s\n', fullfile(out_root,'heat_results.mat'));
end
end

% ================= 工具函数 =================
function v = ensure_col(v), v = v(:); end
function s = sanitize(s), s = regexprep(char(s), '[^a-zA-Z0-9-_]', '_'); end
function f = topfile(ME), if isempty(ME.stack), f='?'; else, f=ME.stack(1).name; end, end
function l = topline(ME), if isempty(ME.stack), l=0; else, l=ME.stack(1).line; end, end

function [variant, cluster] = resolve_labels(Features)
% variant
if ismember('variant',Features.Properties.VariableNames)
    v = string(Features.variant);
else
    v = repmat("Unknown",height(Features),1);
end
variant = v;

% cluster
if ismember('cluster',Features.Properties.VariableNames)
    c = double(Features.cluster);
else
    % 尝试读取同目录聚类文件 byK(K=3)
    c = nan(height(Features),1);
    try
        [fd,~,~] = fileparts(which('counterfactual_main_heat'));
        [in_dir,~,~] = fileparts(which('counterfactual_main_heat')); %#ok<ASGLU>
    catch
    end
    try
        % 假设与 features 同一路径
        [feat_dir,~,~] = fileparts(Features.Properties.Description); %#ok<NASGU>
    catch
    end
    try
        % 最简单：当前工作目录下
        if isfile('smf_fpca_cluster_results.mat')
            Sclu = load('smf_fpca_cluster_results.mat');
            if isfield(Sclu,'byK')
                KK = arrayfun(@(z) z.K, Sclu.byK);
                pos = find(KK==3,1);   % 你说K=3但内部是4原型
                if ~isempty(pos)
                    lbl = Sclu.byK(pos).labels;
                    if numel(lbl)==numel(c)
                        c = double(lbl(:));
                    end
                end
            end
        end
    catch
    end
end
cluster = c;
end

function init = ensure_init(x)
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

function Th = pick_Th_struct(A, t0, T)
post_len = max(0, T - t0 + 1);
if isfield(A,'auto_T') && A.auto_T
    Th = min(A.T_h, max(3, round(A.auto_T_frac * max(1,post_len))));
else
    Th = A.T_h;
end
end

function f_cf = apply_neg_eta(f0, eta)
% 仅放大/缩小负相；eta=-1 => 负相清零；eta< -1 不允许（截到 -1）
eta = max(-1, eta);
scale = 1 + eta;  % ∈ [0, 3]
mneg = (f0 < 0);
f_cf = f0;
f_cf(mneg) = scale * f0(mneg);
end

% === TGC（允许负相，地板相对峰值）===
function f_cf = synthesize_tgc_r2zero_masked_zero(f_base, init0, beta, gamma, prm)
% 目标导向控制（零阈值口径，稳健版）：
%  - 仅在 t>=t0 且“基线为正相”的位置更新 f，使 r(t) 在 T_h 内指数回落；
%  - 允许进入负相，但不低于 -neg_floor，其中 neg_floor = neg_floor_rel * f_peak_pos；
%  - 对基线负相点（t>=t0）可选轻微加深 deepNeg（保持“更负”单调性）；
%  - 全流程严格列向量、边界裁剪，避免维度不一致。

    % ---------- 读参与稳健化 ----------
    f_base = f_base(:);
    T      = numel(f_base);

    t0     = max(1, min(T, round(getfield_def(prm,'t0',1))));
    Th     = max(3, round(getfield_def(prm,'T_h',10)));
    lambda = min(max(getfield_def(prm,'lambda',0.7), 0), 1);
    fbounds= getfield_def(prm,'f_bounds',[-3,3]);
    lo     = fbounds(1); hi = fbounds(2);
    maxIter= max(1, round(getfield_def(prm,'max_iter',3)));

    % 负向地板（相对正峰值的比例；允许进入负相，但不小于 -neg_floor）
    neg_floor_rel = getfield_def(prm,'neg_floor_rel', 0.50); % 0~1
    fpk_pos = max(0, max(f_base));      % 正峰值（若全负则=0）
    neg_floor = neg_floor_rel * max(1e-8, fpk_pos);   % 标量 >=0

    deepNeg = max(0, getfield_def(prm,'deepen_neg', 0.10)); % 仅对基线负相点生效

    % ---------- 初始化 ----------
    f_curr = f_base;
    kdecay = 1 / Th;

    for it = 1:maxIter
        % 1) 用当前 f 积分：得到 r,s,i
        traj = integrate_seir_odeint_style(f_curr, init0, beta, gamma);
        if numel(traj.t) ~= T
            error('tgc_r2zero_zero:early_stop','integrator stopped at t=%d (expect %d)', numel(traj.t), T);
        end
        S = traj.S(:); I = traj.I(:); R = traj.R(:);
        N = S + I + R; N(N<=0) = eps;
        s = S ./ N; i = I ./ N;
        r = beta .* s - gamma;          % 列向量 T×1

        % 2) 构造目标 r_tar（t>=t0 指数回落）
        r_tar = r;                      % 先拷贝
        if t0 <= T
            idx = (t0:T).';             % 列向量
            L   = numel(idx);
            rt0 = r(t0);                % 标量
            % 用 0..L-1 的非负整数生成指数衰减，避免形状不匹配
            dec = exp(-kdecay * ((0:(L-1)).'));   % L×1
            r_tar(idx) = rt0 .* dec;    % 与 idx 等长
        else
            idx = [];                   % 空
        end

        % 3) 反推 f*（全长），随后仅在掩码处更新
        denom  = max(1e-8, 1 - s);      % 避免除 0
        f_star = (N ./ denom) .* ( - (kdecay/beta).*r_tar + beta .* s .* i );
        f_star = f_star(:);

        % 4) 正相更新掩码（基线口径）：仅 t>=t0 且 f_base>0
        mask_update = false(T,1);
        if ~isempty(idx)
            mseg = (f_base(idx) > 0);
            mask_update(idx(mseg)) = true;  % 只打“正相”的点
        end

        f_new = f_curr;

        if any(mask_update)
            seg        = f_curr(mask_update);
            seg_target = f_star(mask_update);

            % 线性混合 + 限幅
            seg = lambda .* seg_target + (1 - lambda) .* seg;
            seg = min(max(seg, lo), hi);

            % 允许负值，但不低于 -neg_floor（固定绝对地板）
            seg = max(seg, -neg_floor);

            % 赋回（两侧长度严格一致）
            if numel(seg) ~= nnz(mask_update)
                error('mask/segment length mismatch: seg=%d, mask=%d', numel(seg), nnz(mask_update));
            end
            f_new(mask_update) = seg;
        end

        % 5) 负相基线：可选轻微加深（仅 t>=t0，且相对“基线”）
        if deepNeg > 0 && ~isempty(idx)
            mneg = (f_base(idx) <= 0);
            if any(mneg)
                idn = idx(mneg);
                f_new(idn) = min(f_new(idn), (1 + deepNeg) .* f_base(idn));
            end
        end

        % 6) 下一轮
        f_curr = f_new;
    end

    f_cf = f_curr(:);
end

function v = getfield_def(S, name, def)
    if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
        v = S.(name);
    else
        v = def;
    end
end



% === 与主文一致的“严格→回退”积分 ===
function traj = safe_integrate(f, init0, beta, gamma, tag)
T = numel(f);
try
    traj = integrate_seir_odeint_style(f, init0, beta, gamma);
    if numel(traj.t) ~= T
        error('safe_integrate:early_stop','t=%d (expect %d)', numel(traj.t), T);
    end
catch
    % 1) 优先用 ode45 放宽
    try
        tspan = 0:(T-1);
        y0    = [init0.Se0; init0.S0; init0.I0; init0.R0];
        opts2 = odeset('RelTol',1e-6,'AbsTol',1e-9);
        [tt, YY] = ode45(@rhs, tspan, y0, opts2);
        if numel(tt) ~= T, error('len'); end
        Se=YY(:,1); S=YY(:,2); I=YY(:,3); R=YY(:,4);
        den = S+I+R; den(den<=0)=eps; I_ode = beta.*S.*I./den;
        traj = struct('t',tt,'Se',Se,'S',S,'I',I,'R',R,'I_ode',I_ode);
        warning('[safe_integrate] relaxed ode45 used for %s', tag);
        return;
    catch
        % 2) 刚性回退 ode15s
        tspan = 0:(T-1);
        y0    = [init0.Se0; init0.S0; init0.I0; init0.R0];
        opts3 = odeset('RelTol',1e-7,'AbsTol',1e-10,'NonNegative',1:4);
        [tt, YY] = ode15s(@rhs, tspan, y0, opts3);
        if numel(tt) ~= T
            error('safe_integrate:still_early_stop','t=%d (expect %d)', numel(tt), T);
        end
        Se=YY(:,1); S=YY(:,2); I=YY(:,3); R=YY(:,4);
        den = S+I+R; den(den<=0)=eps; I_ode = beta.*S.*I./den;
        traj = struct('t',tt,'Se',Se,'S',S,'I',I,'R',R,'I_ode',I_ode);
        warning('[safe_integrate] ode15s fallback used for %s', tag);
        return;
    end
end

    function dydt = rhs(tt, y)
        Se = max(0.0,y(1)); S = max(0.0,y(2)); I = max(0.0,y(3)); R = max(0.0,y(4));
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

function out = integrate_seir_odeint_style(f, init, beta, gamma)
T = numel(f);
tspan = 0:(T-1);
y0 = [max(0,init.Se0); max(0,init.S0); max(0,init.I0); max(0,init.R0)];
    function idx = idx_from_t(tt)
        rr  = round_half_to_even(tt);
        rr  = max(0, min(T-1, rr));
        idx = rr + 1;
    end
    function val = NN_f(tt), val = f(idx_from_t(tt)); end
    function dydt = rhs(tt, y)
        Se = max(0.0,y(1)); S = max(0.0,y(2)); I = max(0.0,y(3)); R = max(0.0,y(4));
        den = S + I + R; if den <= 0, den = eps; end
        inc = beta * S * I / den;
        ft  = NN_f(tt);
        dydt = [-ft; ft - inc; inc - gamma*I; gamma*I];
    end
opts = odeset('RelTol',1e-7,'AbsTol',1e-10,'NonNegative',1:4);
[tt, YY] = ode45(@rhs, tspan, y0, opts);
Se = max(0, YY(:,1)); S = max(0, YY(:,2)); I = max(0, YY(:,3)); R = max(0, YY(:,4));
den = S + I + R; den(den<=0) = eps;
I_ode = beta .* S .* I ./ den;
out = struct('t',tt, 'Se',Se, 'S',S, 'I',I, 'R',R, 'I_ode',I_ode);
end

function r = round_half_to_even(x)
fx=floor(x); frac=x-fx; r=fx;
gt = frac>0.5; eq = abs(frac - 0.5)<=1e-12;
r(gt)=fx(gt)+1; odd = mod(fx,2)~=0; bump = eq & odd; r(bump)=fx(bump)+1;
end
