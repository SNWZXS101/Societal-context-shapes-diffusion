function hbm_uq_plotonly(root_uq, force_rebuild)
% 并列两图：左 \hat{R}（颜色+数字），右 ESS（颜色+数字）
% 
% hbm_uq_plotonly('results_hbm_uq', true);   % 第二个参数=true 强制重建
% 
% hbm_uq_plotonly('results_hbm_uq');
% 




if nargin<1||isempty(root_uq), root_uq='results_hbm_uq'; end
if nargin<2, force_rebuild = false; end

meta_dir    = fullfile(root_uq,'_meta');
cache_stats = fullfile(meta_dir,'uq_stats.mat');
if ~isfolder(meta_dir), mkdir(meta_dir); end

% 统一渲染（稳妥起见仍在轴上单独再设一遍）
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% ==== 读/建统计量（允许强制重建，避免老缓存里的“旧文本”）====
if ~force_rebuild && exist(cache_stats,'file')==2
    S = load(cache_stats);
    Rhat_med   = S.Rhat_med;
    ESS_min    = S.ESS_min;
    tlabels    = S.tlabels;
    spec_names = S.spec_names;
else
    d = dir(root_uq);
    cand = {d([d.isdir]).name};
    cand = cand(~ismember(cand,{'.','..','_meta'}));

    targets = {'delta_early','kappa_cross','rho_eps','tf1_mean_early','tf2_peak','t_tf2_peak'};
    tlabels = lx_targets(targets);
    nT = numel(targets); nS = numel(cand);

    Rhat_med = nan(nT,nS);
    ESS_min  = nan(nT,nS);
    for s=1:nS
        for ti=1:nT
            sumcsv = fullfile(root_uq,cand{s},'M0_main',targets{ti},'summary.csv');
            if exist(sumcsv,'file')~=2, continue; end
            T = readtable(sumcsv);
            if ismember('rhat',T.Properties.VariableNames)
                Rhat_med(ti,s) = median(T.rhat(~isnan(T.rhat)));
            end
            if ismember('ess',T.Properties.VariableNames)
                ESS_min(ti,s)  = min(T.ess(~isnan(T.ess)));
            end
        end
    end
    % ★ 生成“可渲染”的 LaTeX 刻度标签（cellstr）
    spec_names = cellfun(@lx_spec, cand, 'UniformOutput', false);

    save(cache_stats,'Rhat_med','ESS_min','tlabels','spec_names');
end

[nT,nS] = size(Rhat_med);

% ==== 并列两图 ====
fig = figure('Color','w','Position',[80 80 1400 560]);
t = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','loose');

% -- 左：\hat{R}
ax1 = nexttile(t,1);
imagesc(ax1,Rhat_med); axis(ax1,'tight');
colormap(ax1,flipud(parula));
caxis(ax1,[0.98 1.02]);
cb1 = colorbar(ax1);
cb1.Label.String = 'median split-$\hat{R}$';
cb1.Label.Interpreter = 'latex';                  % ★★ 关键：colorbar 标签 LaTeX
cb1.Label.FontSize = 12;

title(ax1,'median split-$\hat{R}$','FontSize',13,'Interpreter','latex');
set(ax1,'XTick',1:nS,'XTickLabel',spec_names,'XTickLabelRotation',35, ...
        'YTick',1:nT,'YTickLabel',tlabels, ...
        'TickLabelInterpreter','latex','FontSize',11,'LineWidth',0.8);
xlim(ax1,[0.5 nS+0.5]); ylim(ax1,[0.5 nT+0.5]);

hold(ax1,'on');
for i=1:nT
    for j=1:nS
        if isfinite(Rhat_med(i,j))
            text(ax1,j,i,sprintf('%.3f',Rhat_med(i,j)), ...
                'HorizontalAlignment','center','VerticalAlignment','middle', ...
                'Color','w','FontWeight','bold','FontSize',11, ...
                'Interpreter','latex');
        end
    end
end
hold(ax1,'off');

% -- 右：ESS
ax2 = nexttile(t,2);
imagesc(ax2,ESS_min); axis(ax2,'tight');
colormap(ax2,parula);
caxis(ax2,[0 max(100, nanmax(ESS_min(:)))]);      % 动态上界
cb2 = colorbar(ax2);
cb2.Label.String = 'min ESS';
cb2.Label.Interpreter = 'latex';                  % ★★
cb2.Label.FontSize = 12;

title(ax2,'min ESS','FontSize',13,'Interpreter','latex');
set(ax2,'XTick',1:nS,'XTickLabel',spec_names,'XTickLabelRotation',35, ...
        'YTick',1:nT,'YTickLabel',tlabels, ...
        'TickLabelInterpreter','latex','FontSize',11,'LineWidth',0.8);
xlim(ax2,[0.5 nS+0.5]); ylim(ax2,[0.5 nT+0.5]);

hold(ax2,'on');
for i=1:nT
    for j=1:nS
        if isfinite(ESS_min(i,j))
            text(ax2,j,i,sprintf('%.0f',ESS_min(i,j)), ...
                'HorizontalAlignment','center','VerticalAlignment','middle', ...
                'Color','w','FontWeight','bold','FontSize',11, ...
                'Interpreter','latex');
        end
    end
end
hold(ax2,'off');

% 总标题
tt = title(t,'Convergence summary: $\hat{R}$ (left) and ESS (right)','FontSize',18,'FontWeight','bold');
set(tt,'Interpreter','latex','FontWeight','bold');                     

end

% ====== LaTeX 标签助手 ======
function L = lx_targets(keys)
L = cell(size(keys));
for i=1:numel(keys)
    switch string(keys{i})                     % ★ 用 keys，而不是 targets
        case "delta_early",     L{i} = '$\Delta_{\mathrm{early}}$';
        case "kappa_cross",     L{i} = '$\kappa_{\mathrm{cross}}$';
        case "rho_eps",         L{i} = '$\rho_{\epsilon}$';
        case "tf1_mean_early",  L{i} = '$f_{1}^{\mathrm{early}}$';
        case "tf2_peak",        L{i} = '$f_{2}^{\mathrm{peak}}$';
        case "t_tf2_peak",      L{i} = '$t_{\mathrm{f2peak}}$';  % 按你要求 t_{f2peak}
        otherwise,              L{i} = strrep(keys{i},'_','\_');
    end
end
end

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
