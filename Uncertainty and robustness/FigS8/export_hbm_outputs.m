function export_hbm_outputs(varargin)
% EXPORT_HBM_OUTPUTS
%  读取 M0_main/<target>/draws_raw.mat 并导出：
%   - summary.csv（含 mean/sd/CI/median + rhat/ess）
%   - ranef_province.csv / ranef_variant.csv
%   - ppc.csv（★ 仍为“变换后尺度”，带 tf 列）
%   - marginal_effects.csv / contrasts.csv
%
%  新增参数：
%   'chains'            : 链数（默认 4）
%   'samples_per_chain' : 每条链保留样本数（默认自动推断 = floor(nRows/chains)）
%   其余参数与旧版一致。

p = inputParser;
addParameter(p,'results_root','results_hbm',@ischar);
addParameter(p,'n_ppc_draws',1000,@(x)isnumeric(x)&&x>0);
addParameter(p,'grid_points',51,@(x)isnumeric(x)&&x>=11);
addParameter(p,'chains',4,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'samples_per_chain',[],@(x)isnumeric(x)&&isscalar(x)&&x>=10);
parse(p,varargin{:});
args = p.Results;

% ---------- 缓存 ----------
cache_path = fullfile(args.results_root,'cache','hbm_cache.mat');
S = load(cache_path); cache = S.cache;
out_root = fullfile(args.results_root,'M0_main');

N = cache.N; P = cache.P; V = cache.V;
fe_names_orig = cache.fe_names_orig;
reg_names = cache.reg_cols; cli_names = cache.cli_cols;

% 连续列数
if isfield(cache,'cont_safe')
    Kc_all = numel(cache.cont_safe);
else
    Kc_all = numel(fe_names_orig) - numel(reg_names) - numel(cli_names);
    if Kc_all < 0, Kc_all = 0; end
end

targets = fieldnames(cache.Y);

for i=1:numel(targets)
    yname = targets{i};
    ydir  = fullfile(out_root, yname);
    if ~isfolder(ydir), warning('[export] Missing dir %s. Skip.', ydir); continue; end

    L = load(fullfile(ydir,'draws_raw.mat'));
    if ~isfield(L,'draws') || isempty(L.draws)
        warning('[export] %s draws empty. Skip.', yname); continue;
    end
    draws = L.draws;
    nRows = height(draws);

    % ------ 推断链参数 ------
    C = args.chains;
    if isempty(args.samples_per_chain)
        Skeep = floor(nRows / C);
    else
        Skeep = args.samples_per_chain;
    end
    nUse = C * Skeep;
    if nUse > nRows, nUse = nRows; Skeep = floor(nUse/C); end
    rowUse = 1:nUse;

    % ------ 工具：汇总 + 诊断 ------
    get_stats = @(x) diag_stats(x(rowUse), C, Skeep);

    % ------ 拆出参数列 ------
    alpha      = pick_simple(draws, {'alpha'});
    sigma      = pick_simple(draws, {'sigma'});
    sigma_prov = pick_simple(draws, {'sigma_prov','tau_prov'});
    sigma_var  = pick_simple(draws, {'sigma_var','tau_var'});
    nu         = pick_simple(draws, {'nu'});

    [beta,  Kuse]  = collect_indexed_fuzzy(draws, 'beta');
    [b_prov,Pgot]  = collect_indexed_fuzzy(draws, 'b_prov');
    [b_var, Vgot]  = collect_indexed_fuzzy(draws, 'b_var');
    [mu,    Ngot]  = collect_indexed_fuzzy(draws, 'mu');

    % 兜底尺寸
    if Pgot < P, tmp = nan(nRows,P); if Pgot>0, tmp(:,1:Pgot)=b_prov; end, b_prov=tmp; end
    if Vgot < V, tmp = nan(nRows,V); if Vgot>0, tmp(:,1:Vgot)=b_var;  end, b_var =tmp; end
    if Ngot < N, tmp = nan(nRows,N); if Ngot>0, tmp(:,1:Ngot)=mu;     end, mu    =tmp; end

    % ---------- 掩码 ----------
    y_tr  = cache.Y.(yname);
    ok    = isfinite(y_tr) & all(isfinite(cache.X_fe),2);
    y_obs = y_tr(ok);
    idx_mu= find(ok);
    N_use = numel(idx_mu);
    mu_use= mu(rowUse, idx_mu);

    % ---------- summary.csv（含 rhat/ess）----------
    rows = {};

    rows{end+1} = make_row('alpha',      get_stats(alpha));
    rows{end+1} = make_row('sigma',      get_stats(sigma));
    rows{end+1} = make_row('sigma_prov', get_stats(sigma_prov));
    rows{end+1} = make_row('sigma_var',  get_stats(sigma_var));
    rows{end+1} = make_row('nu',         get_stats(nu));

    Kname = min(Kuse, numel(fe_names_orig));
    for k=1:Kname
        rows{end+1} = make_row(sprintf('beta:%s', fe_names_orig{k}), get_stats(beta(:,k)));
    end
    for k=Kname+1:Kuse
        rows{end+1} = make_row(sprintf('beta_%d', k), get_stats(beta(:,k)));
    end

    sumtab = vertcat(rows{:});
    writetable(sumtab, fullfile(ydir,'summary.csv'));

    % ---------- ranef_* ----------
    ciP = prctile(b_prov(rowUse,:),[2.5 50 97.5],1);
    ciV = prctile(b_var(rowUse,:), [2.5 50 97.5],1);
    ran_prov = table( ...
        cache.prov_names, ...
        mean(b_prov(rowUse,:),1,'omitnan')', ...
        std(b_prov(rowUse,:),0,1,'omitnan')', ...
        ciP(1,:)', ciP(2,:)', ciP(3,:)', ...
        'VariableNames', {'province','mean','sd','ci2_5','median','ci97_5'});
    ran_var  = table( ...
        string(cache.var_levels(:)), ...
        mean(b_var(rowUse,:),1,'omitnan')', ...
        std(b_var(rowUse,:),0,1,'omitnan')', ...
        ciV(1,:)', ciV(2,:)', ciV(3,:)', ...
        'VariableNames', {'variant','mean','sd','ci2_5','median','ci97_5'});
    writetable(ran_prov, fullfile(ydir,'ranef_province.csv'));
    writetable(ran_var,  fullfile(ydir,'ranef_variant.csv'));

    % ---------- PPC（在“变换后尺度”） ----------
    tf = cache.tf_info.(yname);
    nu_eff = nu; if isempty(nu_eff) || all(~isfinite(nu_eff)), nu_eff = 7 * ones(nUse,1); end
    ssel = randsample(nUse, min(nUse, args.n_ppc_draws));
    yrep_tf = nan(numel(ssel), N_use);
    for s=1:numel(ssel)
        ss = ssel(s);
        tdraw = trnd(nu_eff(ss), [1 N_use]);
        yrep_tf(s,:) = mu_use(ss,:) + sigma(ss) * tdraw;
    end
    pct = prctile(yrep_tf,[10 50 90],1);
    tfcol = repmat(string(tf), N_use, 1);
    ppc = table((1:N_use)', mean(yrep_tf,1)', pct(1,:)', pct(2,:)', pct(3,:)', y_obs(:), tfcol, ...
        'VariableNames', {'row','yhat_mean','yhat_p10','yhat_p50','yhat_p90','y_obs','tf'});
    writetable(ppc, fullfile(ydir,'ppc.csv'));

    % ---------- 边际效应（保持原实现：逆变换） ----------
    gp   = args.grid_points; xgrid = linspace(-2,2,gp);
    Kc   = min(Kc_all, Kuse);
    if Kc >= 1
        me_rows = cell(Kc,1);
        ssel_me = ssel;
        for j=1:Kc
            X0z = zeros(1, Kuse);
            Xs  = repmat(X0z, gp, 1); Xs(:, j) = xgrid(:);
            mu_s = zeros(numel(ssel_me), gp);
            for s=1:numel(ssel_me)
                ss = ssel_me(s);
                mu_s(s,:) = alpha(ss) + Xs * beta(ss,1:Kuse)';
            end
            yrep_me = nan(numel(ssel_me), gp);
            for s=1:numel(ssel_me)
                ss = ssel_me(s);
                tdraw = trnd(nu_eff(ss), [1 gp]);
                yrep_me(s,:) = mu_s(s,:) + sigma(ss)*tdraw;
            end
            yrep_me_raw = local_inverse_transform(yrep_me, tf);
            pctm = prctile(yrep_me_raw,[10 50 90],1);

            if j <= numel(fe_names_orig), namej = string(fe_names_orig{j});
            else,                         namej = "beta_" + j;
            end

            me_rows{j} = table( ...
                repmat(namej, gp,1), xgrid(:), ...
                mean(yrep_me_raw,1)', pctm(1,:)', pctm(2,:)', pctm(3,:)', ...
                'VariableNames', {'covariate','grid_std','yhat_mean','yhat_p10','yhat_p50','yhat_p90'});
        end
        me_tab = vertcat(me_rows{:});
        writetable(me_tab, fullfile(ydir,'marginal_effects.csv'));
    end

    % ---------- 对比 ----------
    start_reg = Kc + 1;
    n_reg = min(numel(reg_names), max(0,Kuse-Kc)); reg_idx = start_reg + (0:n_reg-1);
    start_cli = Kc + n_reg + 1;
    n_cli = min(numel(cli_names), max(0,Kuse-(Kc+n_reg))); cli_idx = start_cli + (0:n_cli-1);

    con_rows = {};
    for r=1:n_reg
        eff = beta(rowUse, reg_idx(r));
        con_rows{end+1} = local_contrast_row(eff, sprintf('region: %s vs East', reg_names{r})); %#ok<AGROW>
    end
    for c=1:n_cli
        eff = beta(rowUse, cli_idx(c));
        con_rows{end+1} = local_contrast_row(eff, sprintf('climate: %s vs Temperate', cli_names{c})); %#ok<AGROW>
    end
    if ~isempty(con_rows)
        writetable(vertcat(con_rows{:}), fullfile(ydir,'contrasts.csv'));
    end
end

fprintf('[export_hbm_outputs] Done.\n');
end

% ====== 工具：把一个参数向量做统计 + 诊断 ======
function s = diag_stats(x, C, S)
% x: nUse×1 向量（按链依次拼接），C: 链数，S: 每链样本
if isempty(x) || all(~isfinite(x))
    s = struct('mean',NaN,'sd',NaN,'ci2_5',NaN,'median',NaN,'ci97_5',NaN,'rhat',NaN,'ess',NaN);
    return;
end
x = x(:);
n = C*S;
if numel(x) < n, n = floor(numel(x)/C)*C; S = n/C; x = x(1:n); end
X = reshape(x, S, C);                % S × C

% 统计量
smean = mean(x,'omitnan');
ssd   = std(x,'omitnan');
ci    = prctile(x,[2.5 50 97.5]);
smean_ch = mean(X,1);                % 1×C
svar_ch  = var(X,0,1);               % 1×C
W = mean(svar_ch);                   % within-chain var
B = S * var(smean_ch, 0, 2);         % between-chain var
var_hat = ((S-1)/S)*W + B/S;
rhat_basic = sqrt(var_hat / W);

% split-Rhat：每链对半
half = floor(S/2);
if half>=10
    Xs = [X(1:half,:) ; X(S-half+1:S,:)];
    Xs = reshape(Xs, half, 2*C);
    smean_ch2 = mean(Xs,1);
    svar_ch2  = var(Xs,0,1);
    W2 = mean(svar_ch2);
    B2 = half * var(smean_ch2,0,2);
    var_hat2 = ((half-1)/half)*W2 + B2/half;
    rhat = sqrt(var_hat2 / W2);
else
    rhat = rhat_basic;
end

% ESS（Geyer 初始正序列的简化实现）
rho_sum = 0; t = 1;
rho_old = 1;
% 去中心化
Xc = X - mean(X,1);
den = sum(sum(Xc.^2));
if den<=0, ess = NaN;
else
    while t <= S-2
        % 自相关（按链求，取平均）
        num = 0; den_ch = 0;
        for c=1:C
            xc = Xc(:,c);
            num = num + sum( xc(1:S-t) .* xc(1+t:S) );
            den_ch = den_ch + sum( xc.^2 );
        end
        rho_t = (num) / (den_ch);
        % 成对截断
        if rho_t < 0, break; end
        rho_sum = rho_sum + 2*rho_t;
        rho_old = rho_t;
        t = t + 1;
        if t>1000, break; end
    end
    n_eff = C*S / (1 + rho_sum);
    ess = max(min(n_eff, C*S), 1);
end

s = struct('mean',smean,'sd',ssd,'ci2_5',ci(1),'median',ci(2),'ci97_5',ci(3), ...
           'rhat',rhat,'ess',ess);
end

function Trow = make_row(name, s)
Trow = table(string(name), s.mean, s.sd, s.ci2_5, s.median, s.ci97_5, s.rhat, s.ess, ...
    'VariableNames', {'param','mean','sd','ci2_5','median','ci97_5','rhat','ess'});
end

% ===== 余下工具（与旧版相同） =====
function x = pick_simple(T, cand)
x = [];
vn = string(T.Properties.VariableNames);
for i=1:numel(cand)
    hit = find(strcmpi(vn, cand{i}), 1); % 不区分大小写
    if ~isempty(hit), x = T.(vn(hit)); return; end
end
end

function [M, count, names] = collect_indexed_fuzzy(T, base)
vars = string(T.Properties.VariableNames);
ord  = nan(size(vars));
keep = false(size(vars));
pat  = "(?i)" + base + ".*?(\d+)(?!.*\d)";
for k=1:numel(vars)
    tok = regexp(vars(k), pat, 'tokens', 'once');
    if ~isempty(tok)
        ord(k) = str2double(tok{1});
        keep(k)= true;
    end
end
names_all = vars(keep);
ord_all   = ord(keep);
[ord_sorted, order] = sort(ord_all);
names = names_all(order);
[~, ia] = unique(ord_sorted, 'stable');
names = names(ia);
ord_sorted = ord_sorted(ia);
count = numel(names);
if count==0, M = []; return; end
Sdraw = height(T);
M = nan(Sdraw, count);
for j=1:count, M(:,j) = T.(names(j)); end
end

function yr = local_inverse_transform(y, tf)
switch lower(tf)
    case 'log',    yr = exp(y);
    case 'abslog', yr = exp(y);
    case 'asinh',  yr = sinh(y);
    otherwise,     yr = y;
end
end

function t = local_contrast_row(eff, label)
ci = prctile(eff,[2.5 50 97.5]);
t = table(string(label), mean(eff,'omitnan'), std(eff,0,'omitnan'), ci(1), ci(2), ci(3), ...
    'VariableNames', {'contrast','mean','sd','ci2_5','median','ci97_5'});
end
