%% ==================== FIG.2b — 原型曲线（中位 ± IQR）直接运行版 · 新版 ====================
% 升级：每面板自适应坐标；4 簇各自颜色；更简洁美观；命令行打印各类样本数（兼容 cluster/byK）。
clear; clc; close all;

%% 0) 用户配置
results_mat  = 'smf_fpca_cluster_results.mat';
outdir       = 'figs_main';

% —— 你期望读取的标签来源 —— 
LABEL_SRC    = 'byK';     % 'cluster' 或 'byK'
K_TARGET     = 3;         % 当 LABEL_SRC='byK' 时使用，例如 3/4
EXPECT_K     = [];        % 可留空[]；或设 4/3 用于强校验

% —— 颜色：优先用你单独的 smf_colorlib；也可在 CLUSTER_COLORS 手动写死 —— 
PALETTE_KEY   = 'classic4';        % 可试：clusters4 / clusters4_alt1/alt2/alt3 / classic4(_deep) 等
CLUSTER_COLORS = [];                     % 若为空，自动从 smf_colorlib(PALETTE_KEY) 取 K 色
% 示例手动：
% CLUSTER_COLORS = [0.10 0.28 0.46; 0.78 0.64 0.00; 0.13 0.55 0.28; 0.75 0.11 0.20];

% —— 样式 —— 
RIBBON_ALPHA = 0.22;     % IQR 带透明
MED_LW       = 2.6;      % 中位线宽
DBA_LW       = 1.6;      % DBA 线宽
DBA_STYLE    = '--';     % DBA 线型
SHOW_DBA     = true;     % 有原型就叠加
SHOW_ZERO    = true;     % 是否画 y=0 细线
GRID_ALPHA   = 0.25;     % 网格浅度
TITLES = {'C1 快起快收','C2 稳起稳收','C3 缓起缓收','C4 多波反弹'};

% —— 导出 —— 
DPI_EXPORT   = 800;

%% 1) 读取并规整（含标签来源解析 + 计数打印）
S = load(results_mat);
[Rbase, labels, K, src_used] = resolve_for_ribbon(S, LABEL_SRC, K_TARGET);
if isempty(EXPECT_K), EXPECT_K = K; end
smf_print_cluster_counts(labels, 'Tag', sprintf('Fig2b:%s', src_used), 'ExpectK', EXPECT_K);

% 基本曲线数据
R = Rbase;   % R.u (M×1), R.F_norm (M×N), R.proto (MxK 或空)

% 输出目录
if ~exist(outdir,'dir'), mkdir(outdir); end

%% 2) 颜色（确保 K 簇各自一色）
C = get_cluster_colors(PALETTE_KEY, K, CLUSTER_COLORS);

%% 3) 作图（每面板自适应坐标）
f = figure('Color','w','Position',[100 100 1200 680]);
tiledlayout(ceil(K/2), 2, 'TileSpacing','compact','Padding','compact');

title_str = sprintf('Prototype ribbons (median ± IQR), K=%d', K);
abox = annotation(f,'textbox',[0 0.965 1 0.035], ...
    'String', title_str, 'EdgeColor','none', 'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', 'FontWeight','bold', 'FontName','Arial', 'FontSize',12);
uistack(abox,'top');

for k = 1:K
    nexttile; hold on; box off; grid on;

    % 柔一点的网格
    ax = gca; ax.GridAlpha = GRID_ALPHA; ax.LineWidth = 0.8;

    idx = (labels==k); n = nnz(idx);
    if n==0
        title(sprintf('Cluster %d (n=0)', k)); continue;
    end

    Fk  = R.F_norm(:, idx);
    med = median(Fk, 2, 'omitnan');
    q25 = quantile(Fk, 0.25, 2);
    q75 = quantile(Fk, 0.75, 2);

    % —— IQR 带 + 中位线 —— 
    fill([R.u; flipud(R.u)], [q25; flipud(q75)], C(k,:), ...
        'FaceAlpha', RIBBON_ALPHA, 'EdgeColor','none');
    plot(R.u, med, '-', 'Color', darken(C(k,:),0.20), 'LineWidth', MED_LW);

    % —— DBA 原型（可选） —— 
    has_proto_k = SHOW_DBA && ~isempty(R.proto) && size(R.proto,2)>=k && size(R.proto,1)==numel(R.u);
    if has_proto_k
        plot(R.u, R.proto(:,k), DBA_STYLE, 'Color', darken(C(k,:),0.35), 'LineWidth', DBA_LW);
    end

    % —— 自适应坐标（每面板） —— 
    if has_proto_k
        ylo = min([q25; R.proto(:,k)]);
        yhi = max([q75; R.proto(:,k)]);
    else
        ylo = min(q25); yhi = max(q75);
    end
    pad = 0.06*(yhi - ylo + eps);
    ylim([ylo - pad, yhi + pad]);

    if SHOW_ZERO
        yl = ylim; plot([R.u(1) R.u(end)], [0 0], '-', 'Color',[0 0 0 0.15], 'LineWidth', 1, 'HandleVisibility','off');
        ylim(yl);  % 保持范围
    end

    xlabel('Normalized time (0–1)'); ylabel('Normalized amplitude');
    if k <= numel(TITLES)
        title(sprintf('%s (n=%d)', TITLES{k}, n));
    else
        title(sprintf('Cluster %d (n=%d)', k, n));
    end
end


% 保存
fname = fullfile(outdir, sprintf('fig2b_prototype_ribbons_K%d', K));
try
    exportgraphics(f, [fname,'.png'], 'Resolution', DPI_EXPORT, 'BackgroundColor','white');
catch
    print(f, [fname,'.png'], '-dpng', ['-r' num2str(DPI_EXPORT)]);
end
savefig(f, [fname,'.fig']);
disp(['Saved: ', fname, '.png 以及 .fig']);

%% ============================ 工具函数 ============================
function [R, labels, K, src_used] = resolve_for_ribbon(S, src, Kreq)
% 解析 u/F_norm/proto 与 labels（支持 cluster/byK）
flds = fieldnames(S);
if numel(flds)==1 && isstruct(S.(flds{1})), S = S.(flds{1}); end
pick = @(T, keys) local_pick(T, keys);

R = struct;
R.u      = pick(S, ["u","grid","time_grid","U"]);
R.F_norm = pick(S, ["F_norm","F","curves","X","Y"]);
if isempty(R.u) || isempty(R.F_norm), error('结果中缺少 u / F_norm'); end
if size(R.F_norm,1) ~= numel(R.u) && size(R.F_norm,2) == numel(R.u)
    R.F_norm = R.F_norm.'; 
end
if isrow(R.u), R.u = R.u.'; end

labels=[]; K=[]; R.proto=[]; src_used='';

switch lower(src)
    case 'byk'
        if isfield(S,'byK') && ~isempty(S.byK)
            KK = arrayfun(@(z) z.K, S.byK);
            pos = find(KK==Kreq,1);
            if ~isempty(pos)
                labels = S.byK(pos).labels;
                K = max(labels);
                R.proto = local_pick(S.byK(pos), ["prototypes","proto","DBA"]);
                src_used = sprintf('byK(K=%d)', Kreq);
            end
        end
    case 'cluster'
        % fallthrough
    otherwise
        % fallthrough
end

if isempty(labels)
    clu = local_struct_pick(S,'cluster');
    labels = pick(clu, ["labels","idx","cluster_id","assignments"]);
    R.proto = pick(clu, ["prototypes","proto","DBA","centers"]);
    if isempty(labels), error('未找到聚类标签 labels'); end
    if isfield(clu,'K') && ~isempty(clu.K), K = clu.K; else, K = max(labels); end
    if isempty(src_used), src_used = 'cluster'; end
end

% 朝向检查：proto需为 MxK
if ~isempty(R.proto)
    if size(R.proto,1)==numel(R.u)
        % ok
    elseif size(R.proto,2)==numel(R.u)
        R.proto = R.proto.'; 
    else
        warning('原型尺寸与 u 不匹配，已跳过 DBA 覆盖。'); R.proto = [];
    end
end
end

function C = get_cluster_colors(palette_key, K, Cmanual)
% 统一拿 K 行颜色
if ~isempty(Cmanual)
    C = Cmanual;
    if size(C,1) < K
        warning('CLUSTER_COLORS 行数不足，已循环补齐。');
        C = repmat(C, ceil(K/size(C,1)), 1);
    end
    C = C(1:K,:); return;
end
% 优先 smf_colorlib
have_lib = exist('smf_colorlib','file')==2;
if have_lib
    try
        C = smf_colorlib(palette_key);
    catch
        C = [];
    end
else
    C = [];
end
% 兜底颜色（四种高区分度）
if isempty(C) || size(C,2)~=3
    C = [0.10 0.28 0.46;   % Navy
         0.78 0.64 0.00;   % Gold
         0.13 0.55 0.28;   % Forest
         0.75 0.11 0.20];  % Crimson
end
% 如颜色少于 K，循环补齐；并仅取前 K 行
if size(C,1) < K
    C = repmat(C, ceil(K/size(C,1)), 1);
end
C = C(1:K,:);
% 若不小心出现重复色，做细微扰动避免“看起来只有两色”的问题
for i=2:K
    if all(abs(C(i,:)-C(i-1,:))<1e-6)
        C(i,:) = min(1, C(i,:).*0.92 + 0.04); % 微调
    end
end
end

function val = local_pick(T, keys)
val = [];
if isempty(T), return; end
if istable(T)
    for i=1:numel(keys)
        k = char(keys(i));
        if any(strcmpi(T.Properties.VariableNames, k)), val = T.(k); return; end
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

function S2 = local_struct_pick(S,k)
S2=[]; if isstruct(S)&&isfield(S,k), S2=S.(k); end
end

function c2 = darken(c, a)
a=max(0,min(1,a)); c2=(1-a)*c;
end

function counts = smf_print_cluster_counts(labels, varargin)
% 打印聚类计数，便于快速核对是否一致（与 S1 / S1cd 版保持一致）
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
