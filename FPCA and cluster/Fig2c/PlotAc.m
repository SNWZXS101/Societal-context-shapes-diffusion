%% ================== FEATURE × PC 横向柱状图（直接运行版·新版） ==================
% 生成两行四块： (c) loadings/associations（PC1, PC2） + (d) contributions%（PC1, PC2）
% - 统一配色：来自 smf_colorlib（含 *_deep 可选），也提供兜底配色
% - 计数打印：兼容 smf_print_cluster_counts（来源：cluster 或 byK(K_TARGET)）
% 输出：figs_ext/s1cd_feature_pc_bars.png 以及 .fig
clear; clc; close all;

%% 0) 用户配置
results_mat  = 'smf_fpca_cluster_results.mat';
features_mat = 'features.mat';
outdir       = 'figs_ext';

TOPK = 5;                         % 显示的特征数量（按 max(|r_PC1|,|r_PC2|) 取 TopK）
FEATURE_WHITELIST = {};           % 如 {'delta_early','rho_eps',...}；空则用全部数值列
MIN_VALID_N = 10;                 % 相关计算最少有效样本（每一列分别判定）

% —— 颜色与样式（来自 colorlib）——
PALETTE_LOAD  = 'velvet_deep';    % (c) 载荷/相关 的颜色主题（建议偏紫等）
PALETTE_CONTR = 'coconut_deep';   % (d) 贡献度 的颜色主题（建议偏青/绿等）
LOAD_IDX      = 3;                % 取主题中的第几个颜色用于 loadings 条
CONTR_IDX     = 1;                % 取主题中的第几个颜色用于 contributions 条
SEP_COLOR     = [0 0 0];          % 0 线颜色
SEP_LW        = 2;                % 0 线宽
BAR_ALPHA     = 0.95;             % 条形不透明度

% —— 计数打印（兼容 smf_print_cluster_counts）——
PRINT_COUNTS = true;              % 是否打印聚类计数
LABEL_SRC    = 'cluster';         % 'cluster' 或 'byK'
K_TARGET     = 4;                 % 当 LABEL_SRC='byK' 时使用（例如 4）
EXPECT_K     = [];                % 可设为空[]，或设为 4/3 等用于核对

% —— 导出 —— 
SAVE_FIG     = true;              % 保存 .fig
SAVE_PNG     = true;              % 保存 .png
DPI_EXPORT   = 600;               % PNG 分辨率
FIG_SIZE_CM  = [22 15];           % 画布（厘米）

%% 1) 读取 scores 与特征表
S = load(results_mat);
[R, Z] = resolve_scores(S);                   % R: meta(ids)；Z: N×PC （取前两列）
[T, X, feat_names, ids_feat] = read_features(features_mat); %#ok<ASGLU>

% —— 顺序对齐（若双方都有 id）——
if ~isempty(R.ids) && ~isempty(ids_feat)
    [ok, ia, ib] = intersect(ids_feat, R.ids, 'stable');
    if isempty(ok)
        warning('features 与 scores 未匹配上 ID，将按原顺序对齐。');
        Zuse = Z; Xuse = X;
    else
        Xuse = X(ia, :);
        Zuse = Z(ib, :);
    end
else
    Zuse = Z; Xuse = X;
end

% —— 可选：只让白名单特征参与 —— 
if ~isempty(FEATURE_WHITELIST)
    [tf, loc] = ismember(string(FEATURE_WHITELIST), feat_names);
    loc = loc(tf);
    if isempty(loc), error('白名单特征未在 features.mat 中找到。'); end
    Xuse = Xuse(:, loc);
    feat_names = feat_names(loc);
end

%% （可选）打印聚类计数，便于每次核对
if PRINT_COUNTS
    [labels_for_count, K_found, src_used] = resolve_labels_for_counts(S, LABEL_SRC, K_TARGET);
    if ~isempty(labels_for_count)
        if isempty(EXPECT_K), EXPECT_K = K_found; end
        smf_print_cluster_counts(labels_for_count, 'Tag',sprintf('S1cd:%s',src_used), 'ExpectK', EXPECT_K);
    else
        warning('S1cd: 未找到可用的 labels（来源=%s）。', LABEL_SRC);
    end
end

%% 2) 与 PC1/PC2 的相关与贡献
[r1, r2, validN] = corr_with_pcs(Zuse, Xuse);    % 行：特征；列：PC1/PC2
% 过滤有效样本不足的特征
keep = validN >= MIN_VALID_N & all(isfinite([r1 r2]),2);
r1 = r1(keep); r2 = r2(keep);
feat_names = feat_names(keep);

% 贡献（按 r^2 归一化到 100%）
contrib1 = 100 * (r1.^2) / max(sum((r1.^2)), eps);
contrib2 = 100 * (r2.^2) / max(sum((r2.^2)), eps);

% 选 TopK（按 max(|r_pc1|, |r_pc2|) 排序）
score_for_rank = max(abs([r1 r2]), [], 2);
[~, ordAll] = sort(score_for_rank, 'descend');
ord = ordAll(1:min(TOPK, numel(ordAll)));

% 子集并倒序（让最强的在最上面）
sel = flipud(ord(:));
names_show = feat_names(sel);
r1_show = r1(sel); r2_show = r2(sel);
c1_show = contrib1(sel); c2_show = contrib2(sel);

%% 3) 颜色（统一来自 colorlib）
C_load_all  = smf_colorlib_safe(PALETTE_LOAD);
C_contr_all = smf_colorlib_safe(PALETTE_CONTR);
LOAD_COLOR  = C_load_all(1+mod(LOAD_IDX-1, size(C_load_all,1)), :);
CONTR_COLOR = C_contr_all(1+mod(CONTR_IDX-1, size(C_contr_all,1)), :);

%% 4) 作图：两行四块（loadings × 2 + contributions × 2）
if ~exist(outdir,'dir'), mkdir(outdir); end
f = figure('Color','w','Units','centimeters','Position',[2 2 FIG_SIZE_CM]);
tiledlayout(f, 2,2,'TileSpacing','compact','Padding','compact');

% ---- (c) PC1 loadings/associations ----
nexttile; hold on; box on; grid on;
barh(r1_show, 'FaceColor', LOAD_COLOR, 'FaceAlpha', BAR_ALPHA, 'EdgeColor','none');
xline(0, '-', 'Color', SEP_COLOR, 'LineWidth', SEP_LW, 'HandleVisibility','off');
set(gca,'YTick',1:numel(names_show),'YTickLabel',pretty_names(names_show), ...
         'YDir','reverse','FontName','Arial','FontSize',18,'LineWidth',1.5, ...
         'TickDir','out','TickLength',[.010 .010]);
xlim([-1 1]); xlabel('Loadings / r'); title('PC1');

% ---- (c) PC2 loadings/associations ----
nexttile; hold on; box on; grid on;
barh(r2_show, 'FaceColor', LOAD_COLOR, 'FaceAlpha', BAR_ALPHA, 'EdgeColor','none');
xline(0, '-', 'Color', SEP_COLOR, 'LineWidth', SEP_LW, 'HandleVisibility','off');
set(gca,'YTick',1:numel(names_show),'YTickLabel',pretty_names(names_show), ...
         'YDir','reverse','FontName','Arial','FontSize',18,'LineWidth',1.5, ...
         'TickDir','out','TickLength',[.010 .010]);
xlim([-1 1]); xlabel('Loadings / r'); title('PC2');

% ---- (d) PC1 contributions ----
nexttile; hold on; box on; grid on;
barh(c1_show, 'FaceColor', CONTR_COLOR, 'FaceAlpha', BAR_ALPHA, 'EdgeColor','none');
set(gca,'YTick',1:numel(names_show),'YTickLabel',pretty_names(names_show), ...
         'YDir','reverse','FontName','Arial','FontSize',18,'LineWidth',1.5, ...
         'TickDir','out','TickLength',[.010 .010]);
xlim([0, max(100, ceil(max(c1_show)*1.15/5)*5)]);
xlabel('Contribution (%)'); title('PC1');

% ---- (d) PC2 contributions ----
nexttile; hold on; box on; grid on;
barh(c2_show, 'FaceColor', CONTR_COLOR, 'FaceAlpha', BAR_ALPHA, 'EdgeColor','none');
set(gca,'YTick',1:numel(names_show),'YTickLabel',pretty_names(names_show), ...
         'YDir','reverse','FontName','Arial','FontSize',18,'LineWidth',1.5, ...
         'TickDir','out','TickLength',[.010 .010]);
xlim([0, max(100, ceil(max(c2_show)*1.15/5)*5)]);
xlabel('Contribution (%)'); title('PC2');

%% 5) 保存
fname = fullfile(outdir,'s1cd_feature_pc_bars');
if SAVE_PNG
    try
        exportgraphics(f, [fname,'.png'], 'Resolution', DPI_EXPORT, 'BackgroundColor','white');
    catch
        print(f, [fname,'.png'], '-dpng', ['-r' num2str(DPI_EXPORT)]);
    end
end
if SAVE_FIG, savefig(f, [fname,'.fig']); end
disp(['Saved: ', fname, '.png', SAVE_FIG*' & .fig']);

%% ============================ 辅助函数 ============================
function [R, Z] = resolve_scores(S)
% 读取前两主成分得分 Z，并返回元信息 R（含 ids）
flds = fieldnames(S);
if numel(flds)==1 && isstruct(S.(flds{1})), S = S.(flds{1}); end
pick = @(T, keys) local_pick(T, keys);
R = struct;
R.ids = local_charcell(pick(S, ["ids","ID","wave_id","waveID","name"]));

Z = pick(S, ["scores","Z","PC_scores","Scores"]);
if isempty(Z), error('结果中未找到 FPCA scores（scores/Z/PC_scores）。'); end
Z = double(Z);
[r,c] = size(Z);
if r < c && r <= 8, Z = Z.'; end  % 修正为 N×PC
if size(Z,2) < 2, error('scores 需要至少两列（PC1/PC2）。'); end
Z = Z(:,1:2);
end

function [T, X, names, ids] = read_features(fmat)
T = load(fmat);
fn = fieldnames(T);
% 取“列最多”的表
Tab = [];
if numel(fn)==1 && istable(T.(fn{1}))
    Tab = T.(fn{1});
else
    best = -inf; picknm = '';
    for i=1:numel(fn)
        if istable(T.(fn{i})) && width(T.(fn{i}))>best
            best = width(T.(fn{i})); picknm = fn{i};
        end
    end
    if ~isempty(picknm), Tab = T.(picknm); end
end
assert(~isempty(Tab), 'features.mat 中未找到 table 变量。');

% id 列（若存在）
ids = string.empty;
idcands = ["wave_id","id","ID","name","waveID"];
hit = intersect(idcands, string(Tab.Properties.VariableNames),'stable');
if ~isempty(hit), ids = string(Tab.(hit(1))); end

% 仅保留数值列
isnum = varfun(@(x) isnumeric(x), Tab, 'OutputFormat','uniform');
TabN = Tab(:, isnum);
names = string(TabN.Properties.VariableNames);
X = table2array(TabN);
end

function [r1, r2, validN] = corr_with_pcs(Z, X)
% 与 PC1/PC2 的 Pearson 相关；对每列单独处理 NaN
Nf = size(X,2);
r1 = nan(Nf,1); r2 = nan(Nf,1); validN = zeros(Nf,1);
x1 = Z(:,1); x2 = Z(:,2);
for j = 1:Nf
    v = X(:,j);
    idx1 = isfinite(v) & isfinite(x1);
    idx2 = isfinite(v) & isfinite(x2);
    if sum(idx1) >= 3
        r1(j) = corr(v(idx1), x1(idx1), 'rows','pairwise', 'Type','Pearson');
    end
    if sum(idx2) >= 3
        r2(j) = corr(v(idx2), x2(idx2), 'rows','pairwise', 'Type','Pearson');
    end
    validN(j) = max(sum(idx1), sum(idx2));
end
% 裁剪到 [-1,1]
r1 = max(-1, min(1, r1));
r2 = max(-1, min(1, r2));
end

function val = local_pick(T, keys)
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

function out = local_charcell(x)
if isempty(x), out = strings(0,1); return; end
if isstring(x), out = x; return; end
if iscellstr(x), out = string(x); return; end
if iscell(x)
    try, out = string(x);
    catch
        out = strings(numel(x),1);
        for i=1:numel(x)
            if ischar(x{i}), out(i)=string(x{i});
            elseif isstring(x{i}), out(i)=x{i};
            else, out(i) = ""; end
        end
    end
elseif ischar(x)
    out = string(x);
else
    out = string(x);
end
end

function nice = pretty_names(s)
% 简化标签：下划线->空格，常见 greek/符号可扩展
nice = regexprep(cellstr(s), '_', ' ');
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
