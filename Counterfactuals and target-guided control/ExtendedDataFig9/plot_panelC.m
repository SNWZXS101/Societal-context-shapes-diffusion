function plot_panelC(out_root, varargin)
% PANEL C — %Avoided 分布（violin / half-violin）+ 关键指标标注
% - 支持按 variant 或 cluster 分组；缺失列会从 Features/Cluster 里自动补齐
% - 单场景：完整小提琴；双场景：左右半提琴
% - 每组标注：n、median、IQR；双场景则两侧分别标注
%
% 依赖：smf_colorlib.m

% ---------- 参数 ----------
p = inputParser; p.FunctionName = mfilename;
p.addRequired('out_root', @(x)ischar(x)||isstring(x));
p.addParameter('GroupBy','variant',@(x)ischar(x)||isstring(x));
p.addParameter('Scenarios',{'NEG_STRONG'},@(x)iscellstr(x)||all(string(x)==string(x)));
p.addParameter('FeaturesMat','',@(x)ischar(x)||isstring(x));
p.addParameter('ClusterMat','',@(x)ischar(x)||isstring(x));
p.addParameter('ClusterKTarget',3,@(x)isnumeric(x)&&isscalar(x));
p.addParameter('GroupOrder',[],@(x)iscellstr(x)||isstring(x)||isempty(x));
p.addParameter('SavePNG',false,@(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('OutName','panelC_violin.png',@(x)ischar(x)||isstring(x));
p.addParameter('FontSize',15,@(x)isnumeric(x)&&isscalar(x)&&x>=6);
p.addParameter('LineWidth',1.2,@(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('FigureSizePx',[1200 700],@(x)isnumeric(x)&&numel(x)==2);
p.addParameter('ShowStats',true,@(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('StatsFontSize',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>0));
p.parse(out_root, varargin{:});
A = p.Results; out_root = char(out_root);
if isempty(A.StatsFontSize), A.StatsFontSize = A.FontSize; end

% ---------- 读取 Results ----------
res_path = fullfile(out_root,'cf_results.mat');
if ~isfile(res_path), error('找不到结果文件：%s（请先运行 counterfactual_main）', res_path); end
S = load(res_path,'Results'); Results = S.Results;

if ~ismember('wave_id', Results.Properties.VariableNames)
    error('Results 中缺少 wave_id，无法与 Features/Clusters 对齐。');
end

% ---------- 若 GroupBy 列缺失则补齐 ----------
groupBy = lower(string(A.GroupBy));
if ~ismember(groupBy, string(Results.Properties.VariableNames))
    switch groupBy
        case "variant"
            if isempty(A.FeaturesMat)
                error('Results 中没有 variant；请通过 ''FeaturesMat'' 提供 features.mat 以补该列。');
            end
            Sf = load(char(A.FeaturesMat),'Features');
            F = Sf.Features(:, {'wave_id','variant'});
            Results = outerjoin(Results, F, 'Keys','wave_id', ...
                                'MergeKeys',true, 'Type','left', ...
                                'LeftVariables', Results.Properties.VariableNames, ...
                                'RightVariables', {'variant'});
        case "cluster"
            if isempty(A.FeaturesMat)
                error('补 cluster 需要 FeaturesMat 以获取行顺序/wave_id。');
            end
            Sf = load(char(A.FeaturesMat),'Features');
            F = Sf.Features(:, {'wave_id'});
            labels = read_cluster_labels(char(A.ClusterMat), A.ClusterKTarget, height(F));
            Tlab = table(string(F.wave_id), labels(:), 'VariableNames', {'wave_id','cluster'});
            Results = outerjoin(Results, Tlab, 'Keys','wave_id', ...
                                'MergeKeys',true, 'Type','left', ...
                                'LeftVariables', Results.Properties.VariableNames, ...
                                'RightVariables', {'cluster'});
        otherwise
            error('未知 GroupBy=%s，须为 variant 或 cluster。', groupBy);
    end
end
if ~ismember(groupBy, string(Results.Properties.VariableNames))
    error('无法补齐分组列 %s。', groupBy);
end

% ---------- 场景 ----------
scns = string(A.Scenarios(:));
if numel(scns)<1 || numel(scns)>2, error('Scenarios 需为 1 或 2 个场景。'); end

% ---------- 组顺序 ----------
groups = unique(string(Results.(groupBy)), 'stable');
if ~isempty(A.GroupOrder)
    go = string(A.GroupOrder(:));
    groups = [go(ismember(go,groups)); groups(~ismember(groups,go))];
end
groups = groups(:).';
groups = ["All" groups];

% ---------- 配色 & 画布 ----------
fs = A.FontSize; sfs = A.StatsFontSize;
f = figure('Color','w','Units','pixels','Position',[120 120 A.FigureSizePx]);
ax = axes('Parent',f); hold(ax,'on'); box(ax,'on'); set(ax,'FontSize',fs);

COLOR = smf_colorlib('waves','K', max(1, numel(scns)));
iqrColor = [0 0 0];
medMarkerSize = 22;

% ---------- 绘制 ----------
X = 1:numel(groups);
xlim(ax,[0.5, numel(groups)+0.5]); grid(ax,'on'); ax.GridAlpha=0.15;
ylabel(ax,'% Avoided'); xlabel(ax, sprintf('Groups by %s', upper(char(groupBy))));
set(ax,'XTick',X,'XTickLabel',groups,'XTickLabelRotation',20);

legend_handles = gobjects(0); legend_texts = strings(0,1);
MAXW = 0.38; OFF = 0.10;                % 半提琴最大半宽，IQR/median 横向偏移
y_text_max = -Inf;                      % 记录标注占用的最大 y，用于扩展 ylim

for g = 1:numel(groups)
    if g==1
        mask_g = true(height(Results),1);
    else
        mask_g = strcmpi(string(Results.(groupBy)), groups(g));
    end

    for s = 1:numel(scns)
        mask_s = strcmpi(string(Results.scenario), scns(s));
        vals = double(Results.avoided_pct(mask_g & mask_s));
        vals = vals(isfinite(vals));
        if isempty(vals), continue; end

        % ---- 小提琴/半提琴 ----
        if numel(scns)==1
            [px, py] = violin_shape(vals, X(g), MAXW);
            ph = patch('XData',px, 'YData',py, 'FaceColor',COLOR, ...
                       'EdgeColor',[0 0 0], 'LineWidth',A.LineWidth, 'FaceAlpha',0.70);
            if g==1, legend_handles(end+1)=ph; legend_texts(end+1)=scns(s); end
        else
            side = 'left'; if s==2, side='right'; end
            [px, py] = half_violin_shape(vals, X(g), MAXW, side);
            ph = patch('XData',px, 'YData',py, 'FaceColor',COLOR(s,:), ...
                       'EdgeColor',[0 0 0], 'LineWidth',A.LineWidth, 'FaceAlpha',0.75);
            if g==1, legend_handles(end+1)=ph; legend_texts(end+1)=scns(s); end
        end

        % ---- 统计：IQR & median ----
        q25 = quantile(vals,0.25); q50 = median(vals); q75 = quantile(vals,0.75);
        if numel(scns)==1
            x_med = X(g);
            plot(ax,[x_med x_med],[q25 q75],'-','Color',iqrColor,'LineWidth',1.4);
            scatter(ax,x_med,q50,medMarkerSize,'MarkerFaceColor','w','MarkerEdgeColor',iqrColor,'LineWidth',1.2);
        else
            x_med = X(g) + (s==1)*(-OFF) + (s==2)*(+OFF);
            plot(ax,[x_med x_med],[q25 q75],'-','Color',iqrColor,'LineWidth',1.2);
            scatter(ax,x_med,q50,medMarkerSize,'MarkerFaceColor','w','MarkerEdgeColor',iqrColor,'LineWidth',1.0);
        end

        % ---- 关键指标标注（n、median、IQR）----
        if A.ShowStats
            n  = numel(vals);
            txt = sprintf('n=%d | med=%.1f%% | IQR=[%.1f, %.1f]%%', n, q50, q25, q75);

            % 单/双场景的文字位置：双场景左右分开、颜色跟随场景色
            if numel(scns)==1
                x_txt = X(g); 
                y_base = max(py);      % 该组小提琴的顶部
                y_txt = y_base + 0.04*(range(py)+eps);
                clr = [0.20 0.20 0.20];
                text(ax, x_txt, y_txt, txt, 'HorizontalAlignment','center','VerticalAlignment','bottom', ...
                    'FontSize', sfs, 'Color', clr, 'Interpreter','none');
            else
                x_txt = X(g) + (s==1)*(-0.22) + (s==2)*(+0.22);
                y_base = max(py);
                y_txt  = y_base + 0.05*(range(py)+eps);
                clr = COLOR(s,:);
                text(ax, x_txt, y_txt, txt, 'HorizontalAlignment','center','VerticalAlignment','bottom', ...
                    'FontSize', sfs, 'Color', clr, 'Interpreter','none');
            end
            y_text_max = max(y_text_max, y_txt);
        end
    end
end

% ---- 标注可能超出坐标，扩展 ylim ----
yl = ylim(ax);
if A.ShowStats && isfinite(y_text_max)
    pad = 0.06*max(1, yl(2)-yl(1));
    ylim(ax, [yl(1), max(yl(2), y_text_max + pad)]);
end

% 图例
if ~isempty(legend_handles)
    lg = legend(ax, legend_handles, legend_texts, 'Location','best', 'Box','off', 'Interpreter','none');
    set(lg,'FontSize',fs);
end

% 标题
if numel(scns)==1
    ttl = sprintf('%%Avoided — %s  |  Grouped by %s', scns(1), upper(char(groupBy)));
else
    ttl = sprintf('%%Avoided — %s vs %s  |  Grouped by %s', scns(1), scns(2), upper(char(groupBy)));
end
title(ax, ttl, 'FontSize', fs+1, 'FontWeight','bold');

% 导出
if A.SavePNG
    fig_dir = fullfile(out_root,'figures'); if ~isfolder(fig_dir), mkdir(fig_dir); end
    out_png = fullfile(fig_dir, char(A.OutName));
    exportgraphics(f, out_png, 'Resolution', 300);
    fprintf('[panelC] Saved → %s\n', out_png);
end
end

% ====== 读取聚类标签（byK优先）======
function labels = read_cluster_labels(cluster_mat, Kreq, Nfeat)
labels = [];
if ~isempty(cluster_mat) && isfile(cluster_mat)
    S = load(cluster_mat);
    flds = fieldnames(S);
    if numel(flds)==1 && isstruct(S.(flds{1})), S = S.(flds{1}); end
    if isfield(S,'byK') && ~isempty(S.byK)
        KK = arrayfun(@(z) getfield_safe(z,'K',NaN), S.byK);
        pos = find(KK==Kreq,1);
        if ~isempty(pos), labels = pickf(S.byK(pos), ["labels","assign","assignments","idx"]); end
    end
    if isempty(labels)
        clu = pickf(S,'cluster');
        labels = pickf(clu, ["labels","assign","assignments","idx"]);
    end
end
if isempty(labels)
    warning('ClusterMat 未提供 cluster 标签；仅当 GroupBy=variant 时可继续。');
    labels = nan(Nfeat,1);
end
labels = labels(:);
end

% ====== 形状函数 ======
function [px, py] = violin_shape(vals, x0, maxw)
y = vals(:); y = y(isfinite(y)); if numel(y)<3, [px,py] = tiny_fallback(y,x0,maxw); return; end
[den, yi] = ksdensity(y, 'Function','pdf');
den = den / max(den) * maxw;
px  = [x0 - den, fliplr(x0 + den)]; py = [yi, fliplr(yi)];
end

function [px, py] = half_violin_shape(vals, x0, maxw, side)
y = vals(:); y = y(isfinite(y)); if numel(y)<3, [px,py] = tiny_fallback(y,x0,maxw); return; end
[den, yi] = ksdensity(y, 'Function','pdf');
den = den / max(den) * maxw;
if strcmpi(side,'left'),  px=[x0, x0 - den, x0]; else, px=[x0, x0 + den, x0]; end
py = [min(yi), yi, max(yi)];
end

function [px,py] = tiny_fallback(y,x0,maxw)
if isempty(y), y = 0; end
lo = min(y); hi = max(y); if lo==hi, lo=lo-0.5; hi=hi+0.5; end
w = maxw*0.15; px = [x0-w x0+w x0+w x0-w]; py = [lo lo hi hi];
end

% ====== 小工具 ======
function val = pickf(T, keys)
val = [];
if isempty(T), return; end
if ischar(keys) || isstring(keys), keys = string(keys); end
if istable(T)
    for i=1:numel(keys)
        k = char(keys(i));
        if any(strcmpi(T.Properties.VariableNames, k)), val = T.(k); return; end
    end
elseif isstruct(T)
        f = fieldnames(T);
        for i=1:numel(keys)
            k = char(keys(i));
            hit = f(strcmpi(f, k)); if ~isempty(hit), val=T.(hit{1}); return; end
        end
end
end

function v = getfield_safe(S, k, def)
v = def; if isstruct(S) && isfield(S,k) && ~isempty(S.(k)), v = S.(k); end
end
