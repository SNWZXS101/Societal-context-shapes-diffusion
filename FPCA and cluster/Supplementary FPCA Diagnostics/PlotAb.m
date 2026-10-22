%% ========================== S1: Scree 直接运行版（含计数打印） ==========================
% 主成分解释率（柱）+ 累计解释率（折线）
% 读取 smf_fpca_cluster_results.mat 中的 expl_var（兼容多种字段名）
% ——新增——
% 读取 labels（cluster / byK），调用 smf_print_cluster_counts 在命令行打印每类样本数
%
% 输出：figs_ext/s1_scree.png 及 .fig
clear; clc; close all;

%% 0) 用户配置
results_mat = 'smf_fpca_cluster_results.mat';
outdir      = 'figs_ext';

% ---- 配色（统一 colorlib）----
PALETTE_KEY = 'classic4_deep';   % 如：clusters4/_alt1/_alt2/_alt3 / waves / sunset / classic4(_deep) 等
BAR_IDX     = 1;                 % 用配色里的第几种色画柱
CUM_IDX     = 4;                 % 用配色里的第几种色画累计线/点

% ---- 样式 ----
MARKER_SIZE = 30;                % 累计点大小
LINE_WIDTH  = 1.8;               % 累计线宽
SHOW_BAR_VALUES = false;         % 在柱顶标注“单个PC解释率”
SHOW_CUM_VALUES = true;          % 在累计点标注“累计解释率”
VALUE_DECIMALS  = 1;             % 百分比小数位
THRESHOLDS  = [80 90];           % 水平虚线阈值（留空 [] 关闭）

% ---- 计数打印（可选）----
PRINT_COUNTS = true;             % 是否打印聚类计数
LABEL_SRC    = 'cluster';        % 'cluster' 或 'byK'
K_TARGET     = 4;                % 当 LABEL_SRC='byK' 时使用（例如 4）
EXPECT_K     = [];               % 可设为空[]，或设为 4/3 等用于核对

% ---- 导出 ----
FIG_SIZE_CM = [18 12];           % 宽×高（厘米）
DPI_EXPORT  = 600;               % PNG 分辨率

%% 1) 读取
S = load(results_mat);

% 1.1 解释率解析
expl = resolve_explained_variance(S);   % 返回列向量（0~1）
assert(~isempty(expl), '未在 %s 中检测到解释率字段（expl_var/explained/vecr 等）。', results_mat);

% 1.2 （可选）打印计数
if PRINT_COUNTS
    [labels_for_count, K_found, src_used] = resolve_labels_for_counts(S, LABEL_SRC, K_TARGET);
    if ~isempty(labels_for_count)
        if isempty(EXPECT_K), EXPECT_K = K_found; end
        smf_print_cluster_counts(labels_for_count, 'Tag',sprintf('S1:%s',src_used), 'ExpectK', EXPECT_K);
    else
        warning('S1: 未找到可用的 labels（来源=%s）。', LABEL_SRC);
    end
end

% 统一为百分比
varpct = 100 * expl(:);
cumpct = cumsum(varpct);
nPC    = numel(varpct);

%% 2) 颜色（统一来自 colorlib）
Call = smf_colorlib_safe(PALETTE_KEY);
BAR_COLOR = Call(1+mod(BAR_IDX-1, size(Call,1)), :);
CUM_COLOR = Call(1+mod(CUM_IDX-1, size(Call,1)), :);

%% 3) 作图
if ~exist(outdir,'dir'), mkdir(outdir); end
f = figure('Color','w','Units','centimeters','Position',[2 2 FIG_SIZE_CM]);
ax = axes(f); hold(ax,'on'); grid(ax,'on'); box(ax,'on');

% 柱
bar(ax, varpct, 'FaceColor', BAR_COLOR, 'EdgeColor', 'none');

% 累计
plot(ax, 1:nPC, cumpct, '-', 'Color', CUM_COLOR, 'LineWidth', LINE_WIDTH);
scatter(ax, 1:nPC, cumpct, MARKER_SIZE, CUM_COLOR, 'filled');

% 阈值线（可选）
ylmax = max([100, ceil(max([cumpct(:); varpct(:)])/5)*5]);
ylim(ax, [0 ylmax]);
if ~isempty(THRESHOLDS)
    for t = THRESHOLDS(:).'
        y = t;
        plot(ax, [0.5 nPC+0.5], [y y], '--', 'Color', [0 0 0 0.35], 'LineWidth', 2, ...
             'HandleVisibility','off');
        text(ax, nPC+0.55, y, sprintf('%d%%', t), 'Color',[0.15 0.15 0.15], ...
             'HorizontalAlignment','left','VerticalAlignment','middle');
    end
end

% 轴与标题
xlabel(ax,'PC');
ylabel(ax,'Explained variance (%)');
title(ax,'FPCA explained variance (bars) and cumulative (line)');
set(ax,'LineWidth',1.5,'FontName','Arial','FontSize',18,'TickDir','out','TickLength',[.010 .010]);

% x 轴刻度
if nPC <= 20
    xticks(ax, 1:nPC);
else
    step = ceil(nPC/10);
    xticks(ax, 1:step:nPC);
end

% 柱顶单PC数值（可选）
if SHOW_BAR_VALUES
    yl = ylim(ax); off = 0.018 * diff(yl);
    for i = 1:nPC
        text(ax, i, varpct(i) + off, sprintf(['%0.' num2str(VALUE_DECIMALS) 'f%%'], varpct(i)), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom', ...
            'FontSize', 18, 'Color',[0.10 0.10 0.10]);
    end
end

% 累计点数值（可选）
if SHOW_CUM_VALUES
    yl = ylim(ax); off = 0.015 * diff(yl);
    for i = 1:nPC
        text(ax, i, cumpct(i) + off, sprintf(['%0.' num2str(VALUE_DECIMALS) 'f%%'], cumpct(i)), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom', ...
            'FontSize', 18, 'Color',CUM_COLOR, 'FontWeight','bold');
    end
end

legend(ax, {'Explained variance','Cumulative'}, 'Location','eastoutside');

%% 4) 高质量导出（PNG + FIG）
fname = fullfile(outdir, 's1_scree');
png_path = [fname,'.png'];
try
    exportgraphics(f, png_path, 'Resolution', DPI_EXPORT, 'BackgroundColor','white');
catch
    print(f, png_path, '-dpng', ['-r' num2str(DPI_EXPORT)]);
end
savefig(f, [fname,'.fig']);
disp(['Saved: ', png_path, ' 以及 ', fname, '.fig']);

%% ============================ 辅助函数 ============================
function expl = resolve_explained_variance(S)
% 在结果结构中鲁棒寻找解释率（返回按比例 0~1 的列向量）
expl = [];
flds = fieldnames(S);
if numel(flds)==1 && isstruct(S.(flds{1})), S = S.(flds{1}); end

% 候选路径：顶层或子结构中常见命名
cands = {{"expl_var"},{"vecr"},{"explained"},{"explained_variance"}};
for i = 1:numel(cands)
    v = local_pick(S, string(cands{i}));
    if ~isempty(v), expl = v; break; end
end

% 若顶层没有，尝试到一层子结构中找
if isempty(expl) && isstruct(S)
    subnames = fieldnames(S);
    for j = 1:numel(subnames)
        T = S.(subnames{j});
        if isstruct(T) || istable(T)
            v = local_pick(T, ["expl_var","vecr","explained","explained_variance"]);
            if ~isempty(v), expl = v; break; end
        end
    end
end

if isempty(expl), return; end
expl = double(expl(:));           % 列向量
if any(expl > 1.001), expl = expl/100; end
expl(expl<0) = 0;
% 归一化（避免轻微数值误差导致累计>100%）
s = sum(expl);
if s > 1.0001, expl = expl/s; end
end

function [labels, K, src_used] = resolve_labels_for_counts(S, src, Kreq)
% 为“打印计数”解析 labels（不影响作图）
labels = []; K = []; src_used = '';
flds = fieldnames(S);
if numel(flds)==1 && isstruct(S.(flds{1})), S = S.(flds{1}); end

switch lower(src)
    case 'cluster'
        if isfield(S,'cluster') && isfield(S.cluster,'labels') && ~isempty(S.cluster.labels)
            labels = S.cluster.labels(:);
            K = max(labels); src_used = 'cluster';
        end
    case 'byk'
        if isfield(S,'byK') && ~isempty(S.byK)
            KK = arrayfun(@(z) z.K, S.byK);
            pos = find(KK==Kreq, 1);
            if ~isempty(pos)
                labels = S.byK(pos).labels(:);
                K = max(labels); src_used = sprintf('byK(K=%d)', Kreq);
            end
        end
    otherwise
        warning('resolve_labels_for_counts: 未知来源 %s（仅支持 cluster/byK）', src);
end
end

function counts = smf_print_cluster_counts(labels, varargin)
% 打印聚类计数，便于快速核对是否一致
p = inputParser;
p.addParameter('Tag','',@(x)ischar(x)||isstring(x));
p.addParameter('ExpectK',[],@(x)isempty(x)||(isscalar(x)&&x==round(x)));
p.parse(varargin{:});
tag = char(p.Results.Tag);  Kexp = p.Results.ExpectK;

labels = labels(:);
bad = ~isfinite(labels);
if any(bad)
    warning('labels 中存在 NaN/Inf，将忽略这些样本计数。');
    labels = labels(~bad);
end
if isempty(labels)
    fprintf('[%s] labels 为空。\n', tag); counts = []; return;
end

K = max(labels);
counts = accumarray(max(1,min(K, round(labels))), 1, [K 1])';
N = numel(labels);

if ~isempty(tag), tinfo = ['[' tag '] ']; else, tinfo = ''; end
fprintf('%sCluster counts (K=%d, N=%d): %s\n', tinfo, K, N, mat2str(counts));
for k = 1:K
    fprintf('  C%d: %d\n', k, counts(k));
end
if ~isempty(Kexp) && K ~= Kexp
    warning('%s实际 K=%d 与期望 K=%d 不一致。', tinfo, K, Kexp);
end
if sum(counts) ~= N
    warning('%s计数之和(%d) 与 N(%d) 不一致。', tinfo, sum(counts), N);
end
end

function val = local_pick(T, keys)
% 在结构体/表 T 中按多个候选名取首个存在的字段/列
val = [];
if isempty(T), return; end
if istable(T)
    for i=1:numel(keys)
        k = char(keys(i));
        if any(strcmpi(T.Properties.VariableNames, k))
            val = T.(k); return;
        end
    end
elseif isstruct(T)
    f = fieldnames(T);
    for i=1:numel(keys)
        k = char(keys(i));
        hit = f(strcmpi(f, k));
        if ~isempty(hit), val = T.(hit{1}); return; end
    end
end
end

function C = smf_colorlib_safe(which)
% 安全调用：若用户未提供 smf_colorlib.m，则给出默认配色
if exist('smf_colorlib','file')==2
    try
        C = smf_colorlib(which);
        return;
    catch
    end
end
% 兜底：Tableau-10
C = [ ...
    0.1216 0.4667 0.7059;
    1.0000 0.4980 0.0549;
    0.1725 0.6275 0.1725;
    0.8392 0.1529 0.1569;
    0.5804 0.4039 0.7412;
    0.5490 0.3373 0.2941;
    0.8902 0.4667 0.7608;
    0.4980 0.4980 0.4980;
    0.7373 0.7412 0.1333;
    0.0902 0.7451 0.8118];
end
