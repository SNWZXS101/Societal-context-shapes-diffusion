function plot_hbm_results(varargin)
% PLOT_HBM_RESULTS
% 可视化导出结果（全英文作图，中文注释）
%
% 图形安排（每个 target 生成若干独立 PNG）：
%  1) Fixed effects — forest (每个 beta 一行；右侧标注 R-hat / ESS)
%  2) Province random intercepts — forest
%  3) Variant random intercepts — forest
%  4) Region & Climate contrasts — two panels
%  5) Posterior predictive check — scatter + 10–90% whiskers
%  6) Marginal effects — 分页（每页最多 10 个 covariates）
%  7) Trace plots — 直接显示 export 生成的 trace PNG（若存在）
%
% 使用：
% plot_hbm_results('results_root','results_hbm', ...
%   'targets',{'delta_early','kappa_cross','rho_eps','xi1','xi2'}, ...
%   'label_map', myLabelMapStruct, ...
%   'dpi',200);
%
% 说明：
% - label_map 用来统一英文名：包含 .covariates / .provinces / .variants / .groups / .targets 五个 Map；
%   不传则使用内置的一个默认 map（可在代码中修改）。

% --------------- 参数 ---------------
p = inputParser;
addParameter(p,'results_root','results_hbm',@ischar);
addParameter(p,'targets',{},@(c)iscell(c)||isstring(c));
addParameter(p,'dpi',200,@(x)isnumeric(x)&&x>50);
addParameter(p,'label_map',[],@(x)isstruct(x)||isempty(x));
parse(p,varargin{:});
args = p.Results;

res_root = args.results_root;
cache_path = fullfile(res_root,'cache','hbm_cache.mat');
S = load(cache_path); cache = S.cache;

% ============ 名称映射（你可在此统一英文名） ============
label_map = args.label_map;
if isempty(label_map)
    label_map = default_label_map();  % 内置一个默认映射
end

% 目标集合
if isempty(args.targets)
    targets = fieldnames(cache.Y);
else
    targets = cellstr(args.targets);
end

out_root = fullfile(res_root,'M0_main');

for ii = 1:numel(targets)
    y = targets{ii};
    ydir = fullfile(out_root, y);
    if ~isfolder(ydir), warning('Missing dir for %s', y); continue; end

    % 读取导出数据
    Sum  = safe_readtable(fullfile(ydir,'summary.csv'));
    PPC  = safe_readtable(fullfile(ydir,'ppc.csv'));
    ME   = safe_readtable(fullfile(ydir,'marginal_effects.csv'));
    RPr  = safe_readtable(fullfile(ydir,'ranef_province.csv'));
    RVa  = safe_readtable(fullfile(ydir,'ranef_variant.csv'));
    Con  = safe_readtable(fullfile(ydir,'contrasts.csv'));
    Diag = safe_readtable(fullfile(ydir,'mcmc_diag.csv'));   % R-hat & ESS

    % FE 名称（支持扩展 FE）
    fe_names = cache.fe_names_orig;
    meta_path = fullfile(ydir,'fe_meta.mat');
    if isfile(meta_path)
        M = load(meta_path);
        if isfield(M,'meta') && isfield(M.meta,'fe_names') && ~isempty(M.meta.fe_names)
            fe_names = M.meta.fe_names;
        end
    end

    % ---------- 1) Fixed effects forest（单张；右侧显示 R-hat/ESS） ----------
    [coef_tab, idx_beta] = get_fixed_effects(Sum, fe_names, label_map);
    % 从诊断表取各 beta 的 Rhat/ESS（列名形式 beta_k）
    beta_diag = [];
    if ~isempty(Diag) && ~isempty(idx_beta)
        % idx_beta 是 beta 的序号向量（与 coef_tab 对齐）
        beta_diag = nan(height(coef_tab), 2);
        for k=1:height(coef_tab)
            j = idx_beta(k);
            row = Diag(strcmp(string(Diag.param), sprintf('beta_%d', j)), :);
            if ~isempty(row)
                beta_diag(k,1) = row.Rhat;
                beta_diag(k,2) = row.ESS;
            end
        end
    end
    if ~isempty(coef_tab)
        f1 = figure('Color','w','Position',[100 100 900 100+24*height(coef_tab)]);
        ax = axes(f1);
        forest_h(ax, coef_tab, 'Effect','95% CI');
        title(sprintf('Fixed effects — %s', label_one(y,label_map,'targets')));
        set(ax,'YDir','reverse'); grid on; xlabel('Posterior mean (with 95% CI)');
        % 在右侧标注 R-hat / ESS
        if ~isempty(beta_diag)
            hold(ax,'on');
            ytick = 1:height(coef_tab);
            xR = max(coef_tab.hi) + 0.05*range([coef_tab.lo;coef_tab.hi]);
            for k=1:numel(ytick)
                if all(isfinite(beta_diag(k,:)))
                    text(ax, xR, ytick(k), sprintf('R̂=%.3f | ESS=%.0f', beta_diag(k,1), beta_diag(k,2)), ...
                        'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',9);
                end
            end
        end
        savepng(f1, fullfile(ydir, sprintf('%s_FE_forest.png', y)), args.dpi);
        close(f1);
    end

    % ---------- 2) Province random intercepts（单张） ----------
    prov_tab = get_ranef_table(RPr, 'province', label_map.provinces);
    if ~isempty(prov_tab)
        f2 = figure('Color','w','Position',[100 100 900 1200]);
        ax2 = axes(f2);
        forest_h(ax2, prov_tab, 'Effect','95% CI');
        title(sprintf('Province random intercepts — %s', label_one(y,label_map,'targets')));
        set(ax2,'YDir','reverse'); grid on; xlabel('Posterior mean (with 95% CI)');
        savepng(f2, fullfile(ydir, sprintf('%s_ranef_province.png', y)), args.dpi);
        close(f2);
    end

    % ---------- 3) Variant random intercepts（单张） ----------
    var_tab = get_ranef_table(RVa, 'variant', label_map.variants);
    if ~isempty(var_tab)
        f3 = figure('Color','w','Position',[100 100 700 450]);
        ax3 = axes(f3);
        forest_h(ax3, var_tab, 'Effect','95% CI');
        title(sprintf('Variant random intercepts — %s', label_one(y,label_map,'targets')));
        set(ax3,'YDir','reverse'); grid on; xlabel('Posterior mean (with 95% CI)');
        savepng(f3, fullfile(ydir, sprintf('%s_variants.png', y)), args.dpi);
        close(f3);
    end

    % ---------- 4) Region & Climate contrasts（并排） ----------
    if ~isempty(Con)
        Creg = Con(startsWith(string(Con.contrast),'region:'),:);
        Ccli = Con(startsWith(string(Con.contrast),'climate:'),:);
        if ~isempty(Creg) || ~isempty(Ccli)
            f4 = figure('Color','w','Position',[100 100 1000 450]);
            tl = tiledlayout(f4,1,2,'TileSpacing','compact','Padding','compact');

            if ~isempty(Creg)
                axL = nexttile(tl,1);
                cont_tab = tidy_contrast(Creg, 'region');
                forest_h(axL, cont_tab, 'Contrast', '95% CI');
                title(axL, sprintf('%s — %s', label_one(y,label_map,'targets'), 'Region (vs East)'));
                set(axL,'YDir','reverse'); grid(axL,'on'); xlabel(axL,'Effect');
            end

            if ~isempty(Ccli)
                axR = nexttile(tl,2);
                cont_tab = tidy_contrast(Ccli, 'climate');
                forest_h(axR, cont_tab, 'Contrast', '95% CI');
                title(axR, sprintf('%s — %s', label_one(y,label_map,'targets'), 'Climate (vs Temperate)'));
                set(axR,'YDir','reverse'); grid(axR,'on'); xlabel(axR,'Effect');
            end

            savepng(f4, fullfile(ydir, sprintf('%s_contrasts_region_climate.png', y)), args.dpi);
            close(f4);
        end
    end

    % ---------- 5) PPC（单张） ----------
    if ~isempty(PPC)
        f5 = figure('Color','w','Position',[100 100 700 500]);
        scatter(PPC.y_obs, PPC.yhat_mean, 20, 'filled'); hold on;
        mn=min(PPC.y_obs); mx=max(PPC.y_obs);
        plot([mn mx],[mn mx],'k--');
        for irow=1:height(PPC)
            plot([PPC.y_obs(irow) PPC.y_obs(irow)], [PPC.yhat_p10(irow) PPC.yhat_p90(irow)], 'Color',[0.6 0.6 0.6]);
        end
        xlabel('Observed'); ylabel('Posterior predictive mean (10–90%)');
        title(sprintf('Posterior predictive check — %s', label_one(y,label_map,'targets'))); grid on;
        savepng(f5, fullfile(ydir, sprintf('%s_ppc.png', y)), args.dpi);
        close(f5);
    end

    % ---------- 6) Marginal effects（分页） ----------
    if ~isempty(ME)
        % 统一协变量英文名（用 label_map）
        ME.covariate = string(ME.covariate);
        for rr=1:height(ME)
            ME.covariate(rr) = label_one(ME.covariate(rr), label_map, 'covariates');
        end
        uCov = unique(string(ME.covariate),'stable');
        per_page = 10; npages = ceil(numel(uCov)/per_page);
        for pg=1:npages
            idx = (pg-1)*per_page + (1:per_page);
            idx = idx(idx<=numel(uCov));
            covs = uCov(idx);

            f6 = figure('Color','w','Position',[100 100 1000 100*max(numel(covs),3)]);
            tl  = tiledlayout(f6, numel(covs), 1, 'TileSpacing','compact','Padding','compact');
            for k=1:numel(covs)
                tab = ME(string(ME.covariate)==covs(k), :);
                ax = nexttile(tl,k);
                patch([tab.grid_std; flipud(tab.grid_std)], [tab.yhat_p10; flipud(tab.yhat_p90)], ...
                      [0.9 0.9 1.0],'EdgeColor','none'); hold on;
                plot(tab.grid_std, tab.yhat_mean, 'b-', 'LineWidth',1.5);
                plot(tab.grid_std, tab.yhat_p50,  'b--');
                xlabel('Std. covariate'); ylabel('Outcome (raw scale)');
                title(char(covs(k)));
                grid on;
            end
            savepng(f6, fullfile(ydir, sprintf('%s_marginals_p%02d.png', y, pg)), args.dpi);
            close(f6);
        end
    end

    % ---------- 7) Trace plots（若已由 export 生成，则直接复制文件名；否则可选跳过） ----------
    % 这里不强画 trace；trace 已在 export 里生成（alpha、sigma、若干 beta）
    % 需要集中排版的话，我也可以改成合并大图；先保留 export 的单图输出即可。

end

fprintf('[plot_hbm_results] Done.\n');

end % main

% ==================== 内部函数 ====================

function T = safe_readtable(fp)
    if isfile(fp), T = readtable(fp); else, T = []; end
end

function s = label_one(name, label_map, group)
% 对单个名进行英文映射（若找不到则保持原样，并把下划线替换为空格）
    if nargin<3, group='covariates'; end
    if isstring(name) || ischar(name), key = string(name);
    else, key = string(name{1}); end
    if isfield(label_map, group) && isa(label_map.(group),'containers.Map') ...
            && isKey(label_map.(group), key)
        s = string(label_map.(group)(key));
    else
        s = replace(key,'_',' ');
    end
end

function [tab, idx_beta] = get_fixed_effects(Sum, fe_names, label_map)
% 从 summary.csv 中抽取 beta 行，并映射英文名；返回顺序与 fe_names 对齐
    tab = []; idx_beta = [];
    if isempty(Sum), return; end
    mask = startsWith(string(Sum.param),'beta:') | startsWith(string(Sum.param),'beta_');
    if ~any(mask), return; end
    T = Sum(mask,:);

    % 将 beta 对应索引与名称配齐
    raw_names = strings(height(T),1);
    idx_beta  = nan(height(T),1);
    for i=1:height(T)
        p = string(T.param(i));
        if startsWith(p,'beta:')
            raw = extractAfter(p,'beta:');         % 原名
            raw_names(i) = raw;
            % 索引：尽量按 fe_names 定位（若能匹配），否则用位置映射
            j = find(strcmp(fe_names, char(raw)), 1);
            if ~isempty(j), idx_beta(i) = j; end
        else
            % 形如 beta_12
            tmp = regexp(p, 'beta_(\d+)$','tokens','once');
            if ~isempty(tmp), idx_beta(i) = str2double(tmp{1}); end
            raw_names(i) = (idx_beta(i) <= numel(fe_names)) * string(fe_names{idx_beta(i)}) ...
                         + (idx_beta(i) >  numel(fe_names)) * string(p);
        end
    end

    labels = cell(height(T),1);
    for i=1:height(T)
        labels{i} = char(label_one(raw_names(i), label_map, 'covariates'));
    end

    tab = table(string(labels), T.mean, T.ci2_5, T.ci97_5, ...
        'VariableNames', {'name','mean','lo','hi'});
end

function tab = get_ranef_table(T, colname, mapobj)
% 读取随机截距表，并做英文名映射
    tab = [];
    if isempty(T) || ~ismember(colname, T.Properties.VariableNames), return; end
    nm = string(T.(colname));
    if isa(mapobj,'containers.Map')
        for i=1:numel(nm)
            k = nm(i);
            if isKey(mapobj,k), nm(i) = mapobj(k); end
        end
    end
    tab = table(nm, T.mean, T.ci2_5, T.ci97_5, ...
       'VariableNames', {'name','mean','lo','hi'});
end

function tab = tidy_contrast(T, which)
% 整理对比表用于作图
    if isempty(T), tab = []; return; end
    nm = string(T.contrast);
    if strcmpi(which,'region')
        nm = replace(nm, "region: ", ""); nm = replace(nm, " vs East", "");
    else
        nm = replace(nm, "climate: ", ""); nm = replace(nm, " vs Temperate", "");
    end
    tab = table(nm, T.mean, T.ci2_5, T.ci97_5, ...
        'VariableNames', {'name','mean','lo','hi'});
end

function forest_h(ax, tab, xlab, slab)
% 通用横向森林图
    axes(ax); cla(ax);
    n = height(tab); if n==0, return; end
    y = 1:n; hold on;
    for i=1:n
        plot([tab.lo(i) tab.hi(i)],[y(i) y(i)],'-','Color',[0.3 0.3 0.3],'LineWidth',2);
    end
    plot(tab.mean, y, 'o','MarkerFaceColor',[0 0.45 0.74],'MarkerEdgeColor','k');
    set(ax,'YTick',y,'YTickLabel',tab.name,'YDir','reverse');
    xlabel(xlab); xline(0,'k:'); box on;
end

function savepng(fig, fp, dpi)
% 安全保存 PNG
    try, exportgraphics(fig, fp, 'Resolution', dpi);
    catch, print(fig, fp, '-dpng', sprintf('-r%d', dpi)); end
end

function L = default_label_map()
% 默认英文映射（你可直接在此修改，或从外部传入 label_map 覆盖）
    L = struct();
    L.covariates = containers.Map( ...
      [ "age65 share","education share","gdp_pc","hospital beds per 1k","household size", ...
        "inbound overnight visitors"," disposable income_{pc}","passenger volume", ...
        "public transit volume","tertiary-sector share","urban population density","urban rate", ...
        "region North","region Northeast","region Central","region South","region Southwest","region Northwest", ...
        "climate Subtropical","climate Tropic","climate Plateau","climate Arid" ...
      ], ...
      [ "Age 65+ share","Education share","GDP per capita","Hospital beds / 1k","Household size", ...
        "Inbound overnight visitors","Disposable income per capita","Passenger volume", ...
        "Public transit volume","Tertiary-sector share","Urban population density","Urbanization rate", ...
        "Region: North","Region: Northeast","Region: Central","Region: South","Region: Southwest","Region: Northwest", ...
        "Climate: Subtropical","Climate: Tropic","Climate: Plateau","Climate: Arid" ...
      ] ...
    );
    L.provinces = containers.Map();  % 如需美化省名，在此添加映射
    L.variants  = containers.Map(["Original","Delta","Omicron"], ["Original","Delta","Omicron"]);
    L.groups    = containers.Map(["region","climate"], ["Region (vs East)","Climate (vs Temperate)"]);
    L.targets   = containers.Map( ...
       ["delta_early","kappa_cross","rho_eps","xi1","xi2"], ...
       ["Early exponential growth","Cross-wave coupling","Reporting noise (asinh)","\xi_1 (global shape)","\xi_2 (global shape)"] );
end
