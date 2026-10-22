function fit_burden_nb_m0_matlab(varargin)
% FIT_BURDEN_NB_M0_MATLAB
% 纯 MATLAB 版 NB2-GLMM（log 链接），两个目标：C_total / I_max
% - 优先使用 反解后的 outcomes: results_hbm/cache/sesir_outcomes.mat
% - 若无，则回退到 features_aug.mat 的 C_total_obs / I_max_obs
% - 支持从上次结果热启动（continue_from_previous=true），或强制从零开始
% - 运行结束自动打印“是否跑好”的体检指标（覆盖率、相关性、ranef幅度、phi范围、热启动参数改变量）
%
% 用法：
% fit_burden_nb_m0_matlab('results_root','results_hbm', ...
%     'features_mat','features_aug.mat', ...
%     'continue_from_previous', true, ...
%     'max_iter', 300, 'max_funevals', 1e5, 'tau_beta', 3.0, ...
%     'n_ppc', 1000, 'display','iter');

p = inputParser;
addParameter(p,'results_root','results_hbm',@(s)ischar(s)||isstring(s));
addParameter(p,'features_mat','features_aug.mat',@(s)ischar(s)||isstring(s));
addParameter(p,'continue_from_previous',true,@islogical);
addParameter(p,'max_iter',300,@(x)isnumeric(x)&&x>0);
addParameter(p,'max_funevals',1e5,@(x)isnumeric(x)&&x>0);
addParameter(p,'tau_beta',3.0,@(x)isnumeric(x)&&x>0);
addParameter(p,'display','iter',@(s)ischar(s)||isstring(s));
addParameter(p,'n_ppc',1000,@(x)isnumeric(x)&&x>0);

% 新增：控制 τ 的下界/是否固定、以及热启动时是否重置 RE
addParameter(p,'tau_lower',0.3,@(x)isnumeric(x)&&x>0);      % τ 的最小值（建议 0.2~0.5）
addParameter(p,'fix_tau',false,@islogical);                 % 是否固定 τp/τv
addParameter(p,'tau_fixed_prov',0.6,@(x)isnumeric(x)&&x>0);
addParameter(p,'tau_fixed_var', 0.6,@(x)isnumeric(x)&&x>0);
addParameter(p,'reinit_random_effects',true,@islogical);    % 热启动时是否重置 b 与 τ

addParameter(p,'offset_mode','none',@(s)ischar(s)||isstring(s));  % 'none'|'logT'|'logCtotal'
addParameter(p,'offset_mode_ctotal','none',@(s)ischar(s)||isstring(s));  % 'none'|'logT'|'logN0'|'logT+logN0'|'logT+logFsum'
addParameter(p,'features_core','features.mat',@(s)ischar(s)||isstring(s)); % 取 init_cell / f_cell 时用


parse(p,varargin{:});
args = p.Results;

% --- 读 HBM 缓存（设计矩阵与索引）
C = load(fullfile(args.results_root,'cache','hbm_cache.mat')); cache = C.cache;
X = cache.X_fe; fe_names = cache.fe_names_orig;
if any(all(abs(X-1)<1e-12,1)), keep = ~all(abs(X-1)<1e-12,1); X = X(:,keep); fe_names = fe_names(keep); end
prov_id = cache.prov_id(:); var_id = cache.var_id(:);
P = cache.P; V = cache.V;

% --- 权重与回退数据：features_aug
Saug = load(args.features_mat);
assert(isfield(Saug,'Features') && istable(Saug.Features), 'features_aug.mat 缺少 Features');
T = Saug.Features;

w = ones(size(X,1),1);
if ismember('w_tail', T.Properties.VariableNames), w = double(T.w_tail);
elseif ismember('w_tail_14', T.Properties.VariableNames), w = double(T.w_tail_14);
end
offset = zeros(size(X,1),1); % 如需按暴露量校正，可改为 log(pop) 或 log(T)

% --- 反解结果（优先）
Ymap = containers.Map;
op = fullfile(args.results_root,'cache','sesir_outcomes.mat');
use_model = isfile(op);
if use_model
    O = load(op); outcomes = O.outcomes;
    if height(outcomes) ~= size(X,1), use_model = false; end
end
if use_model
    Ymap('C_total') = outcomes.C_total_fit(:);
    Ymap('I_max')   = outcomes.I_max_fit(:);
else
    if ismember('C_total_obs', T.Properties.VariableNames), Ymap('C_total') = double(T.C_total_obs(:)); end
    if ismember('I_max_obs',   T.Properties.VariableNames), Ymap('I_max')   = double(T.I_max_obs(:));   end
end

out_root = fullfile(args.results_root,'burden'); if ~isfolder(out_root), mkdir(out_root); end


% ---------- (新增) 计算可选 offset（只对 I_max 使用） ----------
offset_all = zeros(size(X,1),1);
if strcmpi(args.offset_mode,'logT') && ismember('T', T.Properties.VariableNames)
    offset_all = log(max(double(T.T),1));
elseif strcmpi(args.offset_mode,'logCtotal')
    % 需要反解结果
    if exist('outcomes','var')
        offset_all = log(max(outcomes.C_total_fit,1));
    else
        warning('[burden] offset_mode=logCtotal 需要反解 outcomes；未找到 sesir_outcomes.mat，降级为无 offset');
        offset_all = zeros(size(X,1),1);
    end
end

% ---------- (新增) 为 C_total 预计算 offset 组件 ----------
offset_ct_all = zeros(size(X,1),1);

if ~strcmpi(args.offset_mode_ctotal,'none')
    % 取 T（来自 features_aug）
    T_days = zeros(size(X,1),1);
    if ismember('T', T.Properties.VariableNames)
        T_days = double(T.T);
    end

    % 取 N0 与 sum_f（来自 features.mat 的 init_cell / f_cell）
    N0 = zeros(size(X,1),1);
    Fsum = zeros(size(X,1),1);
    try
        Sfeat = load(args.features_core);
        % f_cell
        if isfield(Sfeat,'f_cell'), f_list = Sfeat.f_cell;
        else, f_list = Sfeat.Features.f_cell; end
        % init
        if isfield(Sfeat,'init_cell'), init_list = Sfeat.init_cell;
        else, init_list = Sfeat.Features.init; end

        for n=1:numel(N0)
            ic = init_list{n}; if iscell(ic), ic=ic{1}; end
            Se0 = getfield_fallback(ic, {'Se0','Se','Se_init'}, 0);
            S0  = getfield_fallback(ic, {'S0','S','S_init'},   0);
            I0  = getfield_fallback(ic, {'I0','I','I_init'},   0);
            R0  = getfield_fallback(ic, {'R0','R','R_init'},   0);
            N0(n) = max(S0+I0+R0, 1);                 % N_core 初始
            fsrc  = f_list{n};                         % 外部输入
            Fsum(n) = sum( normalize_f_to_vector(fsrc, max(1,round(T_days(n))) ) );
            % 也可以用 trapz 近似，但离散和足够用了
        end
    catch
        warning('[burden] 无法从 %s 读取 init/f，C_total 的 offset 将退回仅用 T。', args.features_core);
    end

    % 组装 offset_ct_all
    mode = lower(string(args.offset_mode_ctotal));
    if mode=="logt"
        offset_ct_all = log(max(T_days,1));
    elseif mode=="logn0"
        offset_ct_all = log(max(N0,1));
    elseif mode=="logt+logn0"
        offset_ct_all = log(max(T_days,1)) + log(max(N0,1));
    elseif mode=="logt+logfsum"
        offset_ct_all = log(max(T_days,1)) + log(max(Fsum,1));
    else
        offset_ct_all = zeros(size(X,1),1);
    end
end

% —— I_max 的 offset（你之前已加的 logCtotal），保留不动 ——
offset_all_imax = zeros(size(X,1),1);
if strcmpi(args.offset_mode,'logT') && ismember('T', T.Properties.VariableNames)
    offset_all_imax = log(max(double(T.T),1));
elseif strcmpi(args.offset_mode,'logCtotal') && exist('outcomes','var')
    offset_all_imax = log(max(outcomes.C_total_fit,1));
end




for key = ["C_total","I_max"]
    key = char(key);
    if ~isKey(Ymap, key), warning('[burden] 缺少目标 %s，跳过。', key); continue; end
    y_all = Ymap(key);

    % NB2 需要整数：四舍五入并截到非负
    y_all = max(0, round(y_all));

    ok = isfinite(y_all) & (y_all>=0) & (floor(y_all)==y_all) & all(isfinite(X),2);
    if sum(ok) < 20
        warning('[burden] %s 有效样本过少（%d），跳过。', key, sum(ok)); continue;
    end

    y = y_all(ok);
    Xok = X(ok,:);
    wok = max(w(ok),0);

    if strcmp(key,'I_max')
        off = offset_all_imax(ok);
    elseif strcmp(key,'C_total')
        off = offset_ct_all(ok);         % ★ C_total 的 offset
    else
        off = zeros(sum(ok),1);
    end


    pid = prov_id(ok);
    vid = var_id(ok);


    % --- 构造/读取初值（热启动 or 从零）
    ydir = fullfile(out_root, key); if ~isfolder(ydir), mkdir(ydir); end
    init_params = [];
    if args.continue_from_previous
        prev_p = fullfile(ydir,'params.mat');
        if isfile(prev_p)
            P0 = load(prev_p);
            init_params = struct('alpha',P0.alpha,'beta',P0.beta,'bprov',P0.bprov,'bvar',P0.bvar, ...
                'phi',P0.phi,'taup',P0.taup,'tauv',P0.tauv);
            fprintf('[%s] 热启动：使用上次 params.mat 作为初值。\n', key);

            % —— 如果上次是塌陷解（RE≈0 或 τ 很小），或用户要求重置 ——
            if args.reinit_random_effects || (max(abs(P0.bprov))<1e-6 && max(abs(P0.bvar))<1e-6)
                init_params.bprov = 0.05*randn(P,1);
                init_params.bvar  = 0.05*randn(V,1);
                init_params.taup  = max(args.tau_lower, 0.6);
                init_params.tauv  = max(args.tau_lower, 0.6);
                fprintf('[%s] 已重置随机截距与 τ 以跳出塌陷解。\n', key);
            end
        else
            fprintf('[%s] 未发现历史 params.mat，本次从零开始。\n', key);
        end
    else
        fprintf('[%s] 已按要求：从零开始。\n', key);
    end


    % --- 拟合（带 β 岭正则 tau_beta、参数夹断、加大预算）
    [alpha,beta,bprov,bvar,phi,taup,tauv,cov_theta,mu,ll,yrep, fval] = ...
        fit_nb2_map( Xok,y,wok,off,pid,vid,P,V, ...
        args.max_iter,args.max_funevals,args.n_ppc,args.display, ...
        args.tau_beta, args.tau_lower, args.fix_tau, ...
        args.tau_fixed_prov, args.tau_fixed_var, init_params );



    % --- 保存
    save(fullfile(ydir,'params.mat'),'alpha','beta','bprov','bvar','phi','taup','tauv','fe_names','cov_theta');
    save(fullfile(ydir,'ppc.mat'),'mu','ll','y','yrep','pid','vid');

    % --- 体检与打印
    print_diagnostics(key, ydir, y, yrep, bprov, bvar, phi, init_params, alpha,beta,taup,tauv, fval, use_model);
end
end

% ================== 拟合子函数 ==================
function [alpha,beta,bprov,bvar,phi,taup,tauv,cov_theta,mu,ll,yrep,fval] = ...
    fit_nb2_map(X,y,w,offset,pid,vid,P,V,max_iter,max_funevals,n_ppc,disp_flag, ...
    tau_beta, tau_lower, fix_tau, tau_fixed_prov, tau_fixed_var, init_params)


K = size(X,2);
if ~isempty(init_params)
    alpha0 = init_params.alpha; beta0 = init_params.beta(:);
    bprov0 = init_params.bprov(:); bvar0 = init_params.bvar(:);
    phi0 = max(1e-6, init_params.phi); taup0=max(1e-6,init_params.taup); tauv0=max(1e-6,init_params.tauv);
else
    alpha0 = log(max(mean(y),1)); beta0 = zeros(K,1);
    bprov0 = zeros(P,1); bvar0 = zeros(V,1);
    phi0 = 1.0; taup0=0.5; tauv0=0.5;
end

theta0 = [alpha0; beta0; bprov0; bvar0; log(phi0); log(taup0); log(tauv0)];

opts = optimoptions('fminunc','Algorithm','quasi-newton','Display',disp_flag, ...
    'MaxIterations',max_iter,'MaxFunctionEvaluations',max_funevals, ...
    'OptimalityTolerance',1e-6,'StepTolerance',1e-8, ...
    'FiniteDifferenceType','central','SpecifyObjectiveGradient',false);

objfun = @(th) nlp_nb2(th, X, y, w, offset, pid, vid, P, V, ...
    tau_beta, tau_lower, fix_tau, tau_fixed_prov, tau_fixed_var);


[theta_hat, fval, ~, ~, ~, hess] = fminunc(objfun, theta0, opts); %#ok<ASGLU>
[alpha,beta,bprov,bvar,phi,taup,tauv] = unpack(theta_hat,K,P,V);

% 夹断防爆
phi  = min(max(phi ,1e-6),1e6);
if fix_tau
    taup = tau_fixed_prov; tauv = tau_fixed_var;
else
    taup = min(max(taup,tau_lower),1e3);
    tauv = min(max(tauv,tau_lower),1e3);
end

taup = min(max(taup,1e-6),1e3);
tauv = min(max(tauv,1e-6),1e3);

eta = alpha + X*beta + bprov(pid) + bvar(vid) + offset;
mu  = exp(eta);
ll  = w .* nb2_logpmf(y, mu, phi);

Sdraw = min(n_ppc, 2000);
yrep = zeros(Sdraw, numel(y));
for s=1:Sdraw
    rate = gamrnd(phi, mu./phi);
    yrep(s,:) = poissrnd(rate);
end

cov_theta = [];
if ~isempty(hess) && all(isfinite(hess(:))), cov_theta = pinv(hess); end
end

function [nll,grad] = nlp_nb2(th,X,y,w,offset,pid,vid,P,V,tau_beta, tau_lower, fix_tau, tau_fixed_prov, tau_fixed_var)
K = size(X,2);
[alpha,beta,bprov,bvar,phi,taup,tauv] = unpack(th,K,P,V);
phi  = min(max(phi ,1e-6),1e6);
if fix_tau
    taup = tau_fixed_prov; tauv = tau_fixed_var;
else
    taup = min(max(taup,tau_lower),1e3);
    tauv = min(max(tauv,tau_lower),1e3);
end

eta = alpha + X*beta + bprov(pid) + bvar(vid) + offset; mu = exp(eta);
ll  = w .* nb2_logpmf(y, mu, phi);

% 先验/惩罚：若 fix_tau，则不对 τ 求解，只保留 b 与 β 的岭项；否则同原式
if fix_tau
    pen = -0.5*( sum((bprov./taup).^2) + sum((bvar./tauv).^2) + sum((beta./tau_beta).^2) );
else
    pen = -0.5*( sum((bprov./taup).^2) + sum((bvar./tauv).^2) + sum((beta./tau_beta).^2) ) ...
        - P*log(taup) - V*log(tauv);
end
nll = -sum(ll) - pen; grad=[];
end


function v = nb2_logpmf(y, mu, phi)
r = phi; p = r./(r+mu);
p = min(max(p,1e-12), 1-1e-12);
v = gammaln(y+r) - gammaln(r) - gammaln(y+1) + r.*log(p) + y.*log1p(-p);
end

function [alpha,beta,bprov,bvar,phi,taup,tauv] = unpack(theta,K,P,V)
i=1; alpha=theta(i); i=i+1;
beta  = theta(i:i+K-1); i=i+K;
bprov = theta(i:i+P-1); i=i+P;
bvar  = theta(i:i+V-1); i=i+V;
phi   = exp(theta(i)); i=i+1;
taup  = exp(theta(i)); i=i+1;
tauv  = exp(theta(i));
end

% ================== 打印体检 ==================
function print_diagnostics(key, ydir, y, yrep, bprov, bvar, phi, init_params, alpha,beta,taup,tauv, fval, use_model)
yhat_p10 = prctile(yrep,10,1)'; yhat_p50 = prctile(yrep,50,1)'; yhat_p90 = prctile(yrep,90,1)';
yhat_p2  = prctile(yrep,2.5,1)'; yhat_p97 = prctile(yrep,97.5,1)';

cov10_90  = mean(y>=yhat_p10 & y<=yhat_p90);
cov2_5_97 = mean(y>=yhat_p2  & y<=yhat_p97);

% Spearman 相关
rho = corr(y, yhat_p50, 'type','Spearman','rows','pairwise');

% ranef 幅度与 phi 范围
max_bprov = max(abs(bprov)); max_bvar = max(abs(bvar));
phi_ok = (phi>=1e-3 && phi<=1e3);

% 若热启动，测参数改变量
rel_change = NaN;
if ~isempty(init_params)
    v0 = [init_params.alpha; init_params.beta(:); init_params.bprov(:); init_params.bvar(:); ...
        init_params.phi; init_params.taup; init_params.tauv];
    v1 = [alpha; beta(:); bprov(:); bvar(:); phi; taup; tauv];
    rel_change = norm(v1 - v0) / max(norm(v0), 1e-12);
end

% 规则给“PASS/CHECK”
pass_cov = (cov10_90>=0.75) && (cov2_5_97>=0.93);
pass_rho = (rho>=0.7);
pass_re  = (max_bprov<=3 && max_bvar<=3);
pass_all = pass_cov && pass_rho && pass_re && phi_ok;

modeStr = 'obs'; if use_model, modeStr='model'; end
fprintf('\n=== [%s] NB2-GLMM 体检报告（%s） ===\n', key, modeStr);
fprintf('Final objective f(x): %.3f\n', fval);
fprintf('Coverage 10–90%%: %.3f   2.5–97.5%%: %.3f\n', cov10_90, cov2_5_97);
fprintf('Spearman rho(y, yhat_p50): %.3f\n', rho);
fprintf('phi: %.4g   (within [1e-3,1e3]? %s)\n', phi, tf(phi_ok));
fprintf('max|b_prov|: %.3f   max|b_var|: %.3f   (<=3? %s)\n', max_bprov, max_bvar, tf(max_bprov<=3 && max_bvar<=3));
if ~isnan(rel_change)
    fprintf('Warm-start param relative change: %.3e\n', rel_change);
end
fprintf('=== Verdict: %s ===\n\n', ternary(pass_all,'PASS','CHECK'));

% 同步写一个小的 coverage.csv 便于溯源
cov_tab = table(cov10_90, cov2_5_97, rho, phi, max_bprov, max_bvar, ...
    'VariableNames', {'cov10_90','cov2_5_97_5','rho_spearman','phi','max_abs_bprov','max_abs_bvar'});
writetable(cov_tab, fullfile(ydir,'coverage.csv'));

    function s=tf(b), if b, s='OK'; else, s='NO'; end; end
    function s=ternary(cond,a,b), if cond, s=a; else, s=b; end; end
end



% 小工具（别处已存在就不用重复放）
function fvec = normalize_f_to_vector(fsrc, T)
tt = (0:T-1).';
if isnumeric(fsrc)
    f = double(fsrc(:));
    if numel(f)==T, fvec=f; else
        x=linspace(0,T-1,numel(f)).'; fvec=interp1(x,f,tt,'linear','extrap');
    end
elseif isstruct(fsrc)&&isfield(fsrc,'t')&&isfield(fsrc,'f')
    x=double(fsrc.t(:)); y=double(fsrc.f(:)); fvec=interp1(x,y,tt,'linear','extrap');
elseif isa(fsrc,'function_handle')
    try, fvals=arrayfun(@(t) fsrc(t), tt); catch, fvals=arrayfun(@(t) fsrc(0,t), tt); end
    fvec=double(fvals(:));
else, fvec = zeros(T,1);
end
fvec(~isfinite(fvec))=0; fvec=max(fvec,0);
end
function val = getfield_fallback(S, names, defaultVal)
val = defaultVal; if ~isstruct(S), return; end
for i=1:numel(names), if isfield(S,names{i}), val=double(S.(names{i})); return; end, end
end