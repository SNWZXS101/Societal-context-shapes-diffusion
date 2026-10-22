function fig = fig1c_dual_corr_heatmap(results_root, varargin)
% FIG1C_DUAL_CORR_HEATMAP (editable, TeX-ready)
% 上：连续协变量相关；下：目标相关。默认 Spearman/pairwise。
% 新增：统一控制标题/刻度/色标的文本解释器（默认 'tex'）。

p = inputParser;
addParameter(p,'dpi',300,@(x)isnumeric(x)&&x>=72);
addParameter(p,'method','spearman',@(s)any(strcmpi(s,{'spearman','pearson'})));
addParameter(p,'prefer_raw_features',true,@islogical);
addParameter(p,'show_diagonal',false,@islogical);
addParameter(p,'mask_upper',true,@islogical);
addParameter(p,'layout','manual',@(s)any(strcmpi(s,{'manual','tiled'})));
addParameter(p,'label_interpreter','tex',@(s)any(strcmpi(s,{'tex','none','latex'})));
parse(p,varargin{:}); A = p.Results;
if nargin<1||isempty(results_root), results_root='results_hbm'; end
interp = lower(string(A.label_interpreter));

% ===== 数据 =====
S = load(fullfile(results_root,'cache','hbm_cache.mat')); cache = S.cache;
Kc = numel(cache.cont_orig);

cov_names = string(cache.cont_orig(:));
cov_names = truncate_labels(cov_names, 24, interp);

% 连续协变量矩阵
Xc = [];
if A.prefer_raw_features
    feat_mat = fullfile(results_root,'..','features_aug.mat');
    if ~isfile(feat_mat), feat_mat = 'features_aug.mat'; end
    try
        F = load(feat_mat);
        if isfield(F,'Features') && istable(F.Features)
            T = F.Features;
            Xc = nan(height(T), Kc);
            for j=1:Kc
                nm = cache.cont_orig{j};
                if ismember(nm, T.Properties.VariableNames)
                    Xc(:,j) = double(T.(nm));
                end
            end
        end
    catch, Xc=[]; end
end
if isempty(Xc), Xc = cache.X_fe(:,1:Kc); end

Cx = corr_generic(Xc, A.method);

targets = string(fieldnames(cache.Y));
Ymat = nan(cache.N, numel(targets));
for i=1:numel(targets), Ymat(:,i) = double(cache.Y.(targets(i))); end
Ct = corr_generic(Ymat, A.method);

% 聚类排序
ordX = cluster_order(Cx);
ordT = cluster_order(Ct);

% 颜色
cmap = diverge_br();

% ===== 画图 =====
fig = figure('Name','Fig1c Correlations','Color','w','Position',[80 80 1180 960]);

switch lower(A.layout)
case 'manual'
    ax1 = axes('Parent',fig,'Units','normalized','Position',[0.08 0.58 0.84 0.33]); hold(ax1,'on');
    ax2 = axes('Parent',fig,'Units','normalized','Position',[0.08 0.12 0.84 0.33]); hold(ax2,'on');
    safe_pick(ax1); safe_pick(ax2);
otherwise
    tl = tiledlayout(fig,2,1,'Padding','compact','TileSpacing','compact');
    ax1 = nexttile(tl); hold(ax1,'on'); ax2 = nexttile(tl); hold(ax2,'on');
end

% -- 上：协变量 --
plot_corr_panel(ax1, Cx(ordX,ordX), cov_names(ordX), cmap, A, interp);
title(ax1,'Covariate correlations (pairwise rows)','FontSize',14,'Interpreter',char(interp));
safe_pick(get(ax1,'Title'));

% -- 下：目标 --
plot_corr_panel(ax2, Ct(ordT,ordT), targets(ordT), cmap, A, interp);
title(ax2,'Target correlations (transformed scale)','FontSize',13,'Interpreter',char(interp));
safe_pick(get(ax2,'Title'));

% 总标题
ht = sgtitle(fig, 'Correlation structure of inputs and targets','Interpreter',char(interp));
set(ht,'FontSize',16,'FontWeight','bold'); safe_pick(ht);

% 导出
outdir = fullfile(results_root,'M0_main','_figs'); if ~isfolder(outdir), mkdir(outdir); end
if A.dpi>0
    print(fig, fullfile(outdir,'Fig1c_dual_corr_heatmap.png'),'-dpng',sprintf('-r%d',A.dpi));
end

% ===== helpers =====
function plot_corr_panel(ax, C, labels, cmap, A, interp)
    imagesc(ax, C);
    axis(ax,'square'); colormap(ax, cmap); caxis(ax,[-1 1]);

    set(ax,'YDir','reverse');                            % 让第一行在顶部
    % ★ 关键：刻度文本解释器改为 TeX/LaTeX（默认 TeX）
    set(ax,'XTick',1:numel(labels),'XTickLabel',labels,'XTickLabelRotation',45, ...
           'YTick',1:numel(labels),'YTickLabel',labels, ...
           'TickLabelInterpreter',char(interp), ...      % ← 这一行生效
           'Layer','top','Box','off','FontSize',10);

    % 只保留下/左边框
    xl = xlim(ax); yl = ylim(ax);
    line(ax, [xl(1) xl(2)], [yl(2) yl(2)], 'Color',[0 0 0], 'LineWidth',0.6, 'HitTest','off');
    line(ax, [xl(1) xl(1)], [yl(1) yl(2)], 'Color',[0 0 0], 'LineWidth',0.6, 'HitTest','off');

    % 遮上三角（含对角）
    if A.mask_upper
        n = size(C,1);
        for i = 1:n
            for j = i:n
                patch('Parent',ax, ...
                      'XData',[j-0.5 j+0.5 j+0.5 j-0.5], ...
                      'YData',[i-0.5 i-0.5 i+0.5 i+0.5], ...
                      'FaceColor',[1 1 1], 'EdgeColor','none', ...
                      'HitTest','off','PickableParts','none');
            end
        end
    end

    % 对角线（可选）
    if A.show_diagonal
        n = size(C,1);
        plot(ax, 0.5+(1:n), 0.5+(1:n), '-', 'Color',[0 0 0]*0.25, 'LineWidth',0.8, ...
            'HitTest','off','PickableParts','none');
    end

    % 颜色条：也改刻度解释器
    cb = colorbar(ax); cb.Ticks=-1:0.5:1; cb.Box='off';
    set(cb,'TickLabelInterpreter',char(interp));          % ← 颜色条刻度解释器
    safe_pick(cb);
end

function C = corr_generic(M, method)
    switch lower(method)
        case 'spearman', C = corr(M,'Type','Spearman','Rows','pairwise');
        otherwise,       C = corr(M,'Type','Pearson','Rows','pairwise');
    end
    C(isnan(C)) = 0;
end

function ord = cluster_order(C)
    D = 1 - C; D = (D + D')/2; D(1:size(D,1)+1:end)=0;
    Y = squareform(max(D,0));
    Z = linkage(Y,'average');
    ord = optimalleaforder(Z, Y);
end

function cmap = diverge_br()
    neg=[0.10 0.35 0.70]; neu=[0.92 0.92 0.92]; pos=[0.82 0.10 0.12];
    t = linspace(0,1,256)'; cmap=zeros(256,3);
    for k=1:256
        s=2*t(k)-1;
        if s>=0, cmap(k,:)=neu+(pos-neu)*s; else, cmap(k,:)=neu+(neg-neu)*(-s); end
    end
end

function s = truncate_labels(a,maxch,interp)
    s=strings(size(a));
    for i=1:numel(a)
        t=char(a(i));
        if numel(t)>maxch, t=[t(1:maxch-1) '...']; end % 用 '...' 兼容 TeX
        % TeX 模式下，如果你不想 '_' 被当成下标，可以在这里转义：
        if strcmpi(interp,'tex')
            t = regexprep(t,'_','\_');
        end
        s(i)=string(t);
    end
end

function safe_pick(h)
    try
        if isgraphics(h)
            if isprop(h,'HitTest'), set(h,'HitTest','on'); end
            if isprop(h,'PickableParts'), set(h,'PickableParts','all'); end
        end
    catch, end
end

end
