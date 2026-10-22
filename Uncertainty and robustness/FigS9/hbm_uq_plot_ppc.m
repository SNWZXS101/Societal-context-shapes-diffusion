function hbm_uq_plot_ppc(root_uq, force_rebuild)
% 并列两图：左 ρ_Spearman（PPC 排名相关），右 10–90% 覆盖率
% 读取 <root_uq>/<spec>/M0_main/<target>/ppc.csv
% 
% hbm_uq_plot_ppc('results_hbm_uq', true);
% 
% 
% hbm_uq_plot_ppc('results_hbm_uq');


if nargin<1||isempty(root_uq), root_uq='results_hbm_uq'; end
if nargin<2, force_rebuild=false; end

meta_dir = fullfile(root_uq,'_meta'); if ~isfolder(meta_dir), mkdir(meta_dir); end
cache_fp = fullfile(meta_dir,'ppc_stats.mat');

% 统一渲染
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

targets = {'delta_early','kappa_cross','rho_eps','tf1_mean_early','tf2_peak','t_tf2_peak'};
tlabels = lx_targets(targets);   % LaTeX y 轴标签

if ~force_rebuild && exist(cache_fp,'file')==2
    S = load(cache_fp);
    RHO = S.RHO; COV = S.COV; spec_names = S.spec_names;
else
    % 列出所有扰动规格
    d = dir(root_uq);
    cand = {d([d.isdir]).name};
    cand = cand(~ismember(cand,{'.','..','_meta'}));
    nT = numel(targets); nS = numel(cand);
    RHO = nan(nT,nS);
    COV = nan(nT,nS);

    for s = 1:nS
        for ti = 1:nT
            ppc_csv = fullfile(root_uq,cand{s},'M0_main',targets{ti},'ppc.csv');
            if exist(ppc_csv,'file')~=2, continue; end
            T = readtable(ppc_csv);
            % Spearman rank correlation（yhat_mean vs y_obs）
            yhat = T.yhat_mean; yobs = T.y_obs;
            msk = isfinite(yhat) & isfinite(yobs);
            if nnz(msk)>=3
                RHO(ti,s) = corr(yhat(msk), yobs(msk), 'Type','Spearman');
            end
            % 10–90% 覆盖率
            lo = T.yhat_p10; hi = T.yhat_p90;
            m2 = isfinite(lo) & isfinite(hi) & isfinite(yobs);
            if nnz(m2)>0
                COV(ti,s) = mean(yobs(m2) >= lo(m2) & yobs(m2) <= hi(m2));
            end
        end
    end
    % 横轴标签（LaTeX）
    spec_names = cellfun(@lx_spec, cand, 'uni', 0);
    save(cache_fp,'RHO','COV','spec_names');
end

[nT,nS] = size(RHO);

% ===== 并列两图 =====
fig = figure('Color','w','Position',[80 80 1400 560]);
t = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','loose');

% 左：Spearman ρ
ax1 = nexttile(t,1);
imagesc(ax1,RHO); axis(ax1,'tight');
colormap(ax1,parula);
caxis(ax1,[-0.1 1.0]);
cb1 = colorbar(ax1);
cb1.Label.String = '$\rho_{\mathrm{Spearman}}$';
cb1.Label.Interpreter = 'latex';
cb1.Label.FontSize = 14;
title(ax1,'Spearman rank correlation','Interpreter','latex','FontSize',14);
set(ax1,'XTick',1:nS,'XTickLabel',spec_names,'XTickLabelRotation',35, ...
        'YTick',1:nT,'YTickLabel',tlabels,'TickLabelInterpreter','latex', ...
        'FontSize',11,'LineWidth',0.8);
xlim(ax1,[0.5 nS+0.5]); ylim(ax1,[0.5 nT+0.5]);
hold(ax1,'on');
for i=1:nT
    for j=1:nS
        if isfinite(RHO(i,j))
            text(ax1,j,i,sprintf('%.2f',RHO(i,j)), ...
                'HorizontalAlignment','center','VerticalAlignment','middle', ...
                'Color','w','FontWeight','bold','FontSize',14,'Interpreter','latex');
        end
    end
end
hold(ax1,'off');

% 右：10–90% 覆盖率
ax2 = nexttile(t,2);
imagesc(ax2,COV); axis(ax2,'tight');
colormap(ax2,parula);
caxis(ax2,[0.6 1.0]);   % 正常应接近 0.8
cb2 = colorbar(ax2);
cb2.Label.String = '10--90\% coverage';
cb2.Label.Interpreter = 'latex';
cb2.Label.FontSize = 14;
title(ax2,'Posterior predictive coverage (10--90\%)','Interpreter','latex','FontSize',13);
set(ax2,'XTick',1:nS,'XTickLabel',spec_names,'XTickLabelRotation',35, ...
        'YTick',1:nT,'YTickLabel',tlabels,'TickLabelInterpreter','latex', ...
        'FontSize',11,'LineWidth',0.8);
xlim(ax2,[0.5 nS+0.5]); ylim(ax2,[0.5 nT+0.5]);
hold(ax2,'on');
for i=1:nT
    for j=1:nS
        if isfinite(COV(i,j))
            text(ax2,j,i,sprintf('%.2f',COV(i,j)), ...
                'HorizontalAlignment','center','VerticalAlignment','middle', ...
                'Color','w','FontWeight','bold','FontSize',14,'Interpreter','latex');
        end
    end
end
hold(ax2,'off');

% 总标题
tt = title(t,'Posterior predictive checks: $\rho$ (left) and coverage (right)','FontSize',18);
set(tt,'Interpreter','latex');

end

% ===== 辅助：LaTeX 目标标签 =====
function L = lx_targets(keys)
L = cell(size(keys));
for i=1:numel(keys)
    switch string(keys{i})
        case "delta_early",     L{i} = '$\Delta_{\mathrm{early}}$';
        case "kappa_cross",     L{i} = '$\kappa_{\mathrm{cross}}$';
        case "rho_eps",         L{i} = '$\rho_{\epsilon}$';
        case "tf1_mean_early",  L{i} = '$f_{1}^{\mathrm{early}}$';
        case "tf2_peak",        L{i} = '$f_{2}^{\mathrm{peak}}$';
        case "t_tf2_peak",      L{i} = '$t_{\mathrm{f2peak}}$';
        otherwise,              L{i} = strrep(keys{i},'_','\_');
    end
end
end

% ===== 辅助：LaTeX 规格刻度 =====
function s = lx_spec(name_raw)
switch string(name_raw)
    case "base",        s = 'base';
    case "no_HS",       s = 'no $HS$';
    case "no_RS",       s = 'no $RS$';
    case "no_RS_HS",    s = 'no $RS{+}HS$';
    case "priorsx0p5",  s = '$\mathrm{priors}\!\times\!0.5$';
    case "priorsx2",    s = '$\mathrm{priors}\!\times\!2$';
    otherwise,          s = strrep(name_raw,'_','\_');
end
end
