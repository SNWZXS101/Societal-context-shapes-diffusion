function export_burden_outputs_matlab(varargin)
% EXPORT_BURDEN_OUTPUTS_MATLAB
% 导出负担层结果（两目标各一套）：
%   - summary.csv：alpha/beta/phi/随机截距
%   - ranef_province.csv / ranef_variant.csv
%   - ppc.csv：预测均值与分位数（基于 NB2 的后验样本模拟）
%
% 兼容目录命名：
%   - 新命名：results_hbm/burden/C_total/ 与 /I_max/
%   - 旧命名：results_hbm/burden/C_total_obs/ 与 /I_max_obs/

p = inputParser;
addParameter(p,'results_root','results_hbm',@(s)ischar(s)||isstring(s));
addParameter(p,'n_ppc_draws',1000,@(x)isnumeric(x)&&x>0);
parse(p,varargin{:});
args = p.Results;

% 读取 FE/索引缓存（取特征名）
cache_path = fullfile(args.results_root,'cache','hbm_cache.mat');
C = load(cache_path); cache = C.cache;

targets_new = {'C_total','I_max'};
targets_old = {'C_total_obs','I_max_obs'};
out_root = fullfile(args.results_root,'burden');

for i=1:numel(targets_new)
    yname_new = targets_new{i};
    yname_old = targets_old{i};

    ydir = fullfile(out_root, yname_new);
    if ~isfolder(ydir)
        % 回退到旧命名
        ydir_old = fullfile(out_root, yname_old);
        if isfolder(ydir_old)
            ydir = ydir_old;
        else
            warning('[export-burden] 缺少目录 %s 或 %s，跳过。', ...
                    fullfile(out_root,yname_new), fullfile(out_root,yname_old));
            continue;
        end
    end

    % 读取抽样（或MAP模拟）输出
    p_params = fullfile(ydir,'params.mat');
    p_ppc    = fullfile(ydir,'ppc.mat');
    if ~isfile(p_params) || ~isfile(p_ppc)
        warning('[export-burden] %s 缺少 params.mat 或 ppc.mat，跳过。', ydir);
        continue;
    end
    P = load(p_params);
    Q = load(p_ppc);

    % ---- summary.csv（alpha/beta/phi/tau*）----
    alpha = P.alpha; beta = P.beta(:); phi = P.phi;
    taup  = P.taup;  tauv = P.tauv;
    fe_names = string(P.fe_names(:));

    sumtab = table();
    sumtab = [sumtab; one_row('alpha', alpha)];
    sumtab = [sumtab; one_row('phi',   phi  )];
    sumtab = [sumtab; one_row('taup',  taup )];
    sumtab = [sumtab; one_row('tauv',  tauv )];
    for k=1:numel(beta)
        sumtab = [sumtab; one_row("beta:"+fe_names(k), beta(k))]; %#ok<AGROW>
    end
    writetable(sumtab, fullfile(ydir,'summary.csv'));

    % ---- ranef ----------
    bprov = P.bprov(:); bvar = P.bvar(:);
    ran_prov = table((1:numel(bprov))', bprov, ...
        'VariableNames', {'province_id','effect'});
    ran_var  = table({'Original';'Delta';'Omicron'}, bvar, ...
        'VariableNames', {'variant','effect'});
    writetable(ran_prov, fullfile(ydir,'ranef_province.csv'));
    writetable(ran_var,  fullfile(ydir,'ranef_variant.csv'));

    % ---- ppc.csv（已在拟合阶段生成 yrep；这里只做汇总）----
    y_obs = Q.y(:);
    yrep  = Q.yrep;           % S x N
    yhat_mean = mean(yrep,1)'; 
    q = prctile(yrep,[10 50 90],1)';   % N x 3
    ppc = table((1:numel(y_obs))', double(y_obs), yhat_mean, q(:,1), q(:,2), q(:,3), ...
        'VariableNames', {'row','y_obs','yhat_mean','yhat_p10','yhat_p50','yhat_p90'});
    writetable(ppc, fullfile(ydir,'ppc.csv'));

    fprintf('[export_burden_outputs_matlab] 导出完成：%s\n', ydir);
end
end

% ===== helpers =====
function tab = one_row(name, val)
tab = table(string(name), double(val), NaN, NaN, NaN, NaN, ...
    'VariableNames', {'param','mean','sd','ci2_5','median','ci97_5'});
end
