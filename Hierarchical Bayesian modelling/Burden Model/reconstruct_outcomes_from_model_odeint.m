function outcomes = reconstruct_outcomes_from_model_odeint(varargin)
% 使用 integrate_seir_odeint_style 逐波反解 SeSIR 轨迹，生成：
%   C_total_fit = ∫ incidence(t) dt = ∫ I_ode(t) dt
%   I_max_fit   = max_t incidence(t) = max_t I_ode(t)
% 说明：
%   - out.I_ode 在 ODE 里定义为 beta * S * I / (S+I+R)，即“新增曲线”（发生率）
%   - 不使用状态量 I(t) 作为峰值或积分

p = inputParser;
addParameter(p,'features_mat','features.mat',@(s)ischar(s)||isstring(s));
addParameter(p,'features_aug','features_aug.mat',@(s)ischar(s)||isstring(s));
addParameter(p,'out_dir','results_hbm',@(s)ischar(s)||isstring(s));
parse(p,varargin{:});
args = p.Results;

S1 = load(args.features_mat);
assert(isfield(S1,'Features') && istable(S1.Features), 'features.mat 缺少 table 变量 Features');
Tfeat = S1.Features;

S2 = load(args.features_aug);
assert(isfield(S2,'Features') && istable(S2.Features), 'features_aug.mat 缺少 table 变量 Features');
Taug = S2.Features;

% ---------- 兼容顶层/表内的 f & init ----------
% f_cell
if isfield(S1,'f_cell')
    f_list = S1.f_cell;   % 210x1 cell
else
    assert(ismember('f_cell', Tfeat.Properties.VariableNames), 'features 缺少 f_cell（顶层或表列）');
    f_list = Tfeat.f_cell;
end
% init
if isfield(S1,'init_cell')
    init_list = S1.init_cell;   % 210x1 cell
else
    assert(ismember('init', Tfeat.Properties.VariableNames), 'features 缺少 init（顶层 init_cell 或表列 init）');
    init_list = Tfeat.init;
end

% 行数一致性
assert(numel(f_list)==height(Tfeat) && numel(init_list)==height(Tfeat), ...
    'f_cell/init_cell 与 Features 行数不一致');

% 与 aug 对齐
assert(height(Tfeat)==height(Taug), 'features 与 features_aug 行数不一致');
N = height(Tfeat);

% 检查 aug 必要列
for nm = ["beta","gamma","T"]
    assert(ismember(nm, Taug.Properties.VariableNames), 'features_aug 缺少列：%s', nm);
end

C_total_fit = zeros(N,1);
I_max_fit   = zeros(N,1);
I_curve     = cell(N,1);   % 状态 I(t)（保留以便需要）
inc_curve   = cell(N,1);   % 新增曲线 I_ode(t)

for n = 1:N
    T = max(1, round(double(Taug.T(n))));
    beta  = double(Taug.beta(n));
    gamma = double(Taug.gamma(n));

    % f：统一为长度 T 的列向量
    fsrc = f_list{n};
    fvec = normalize_f_to_vector(fsrc, T);

    % init：允许是 cell 中的 struct
    ic = init_list{n};
    if iscell(ic), ic = ic{1}; end
    init.Se0 = getfield_fallback(ic, {'Se0','Se','Se_init'}, 0);
    init.S0  = getfield_fallback(ic, {'S0','S','S_init'},   0);
    init.I0  = getfield_fallback(ic, {'I0','I','I_init'},   0);
    init.R0  = getfield_fallback(ic, {'R0','R','R_init'},   0);

    out = integrate_seir_odeint_style(fvec, init, beta, gamma);

    % === 指标定义（基于新增曲线 I_ode）===
    C_total_fit(n) = trapz(out.t, out.I_ode);  % 累计新增（积分）
    I_max_fit(n)   = max(out.I_ode);           % 新增峰值（max incidence）

    I_curve{n}   = out.I(:);
    inc_curve{n} = out.I_ode(:);
end

outcomes = table(C_total_fit, I_max_fit, inc_curve, I_curve);

out_dir = fullfile(args.out_dir,'cache');
if ~isfolder(out_dir), mkdir(out_dir); end
save(fullfile(out_dir,'sesir_outcomes.mat'), 'outcomes');

fprintf('[reconstruct_outcomes_from_model_odeint] OK. N=%d（I_max=peak of incidence, C_total=integral of incidence）\n', N);
end

% ---------- helpers ----------
function fvec = normalize_f_to_vector(fsrc, T)
    tt = (0:T-1).';
    if isnumeric(fsrc)
        f = double(fsrc(:));
        if numel(f)==T
            fvec = f;
        else
            x = linspace(0, T-1, numel(f)).';
            % linspace 生成唯一采样点，避免 interp1 的“采样点必须唯一”报错
            fvec = interp1(x, f, tt, 'linear', 'extrap');
        end
    elseif isstruct(fsrc) && isfield(fsrc,'t') && isfield(fsrc,'f')
        x = double(fsrc.t(:)); y = double(fsrc.f(:));
        assert(numel(x)==numel(y), 'f_cell 结构体 t 与 f 长度不一致');
        % 若存在重复 t，可做轻微抖动避免唯一性报错
        [x, ia] = unique(x,'stable'); y = y(ia);
        fvec = interp1(x, y, tt, 'linear', 'extrap');
    elseif isa(fsrc,'function_handle')
        try, fvals = arrayfun(@(t) fsrc(t), tt);
        catch, fvals = arrayfun(@(t) fsrc(0,t), tt);
        end
        fvec = double(fvals(:));
    else
        error('不支持的 f_cell 类型（向量/struct/function_handle）');
    end
    fvec(~isfinite(fvec)) = 0;
    fvec = max(fvec, 0);
end

function val = getfield_fallback(S, names, defaultVal)
    val = defaultVal;
    if ~isstruct(S), return; end
    for i=1:numel(names)
        if isfield(S, names{i})
            val = double(S.(names{i}));
            return;
        end
    end
end
