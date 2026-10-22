function out = FIG2d_FromMainDTW()
% ...（前半段同你当前版本，不变）...

%% 路径
results_mat = 'smf_fpca_cluster_results.mat';
outdir = 'figs_main'; if ~exist(outdir,'dir'), mkdir(outdir); end

%% 读取 & 规整
S = load(results_mat);
[labels, K, D] = get_labels_and_D(S);   % K= max(labels)
smf_print_cluster_counts(labels, 'Tag','Fig2d', 'ExpectK', K);

%% 计算分离度（与主脚本(d)一致）
intra = nan(K,1);
inter_near = nan(K,1);
nearestCluster = nan(K,1);
for k = 1:K
    Ik = find(labels==k);
    if numel(Ik) >= 2
        Dk = D(Ik, Ik);
        intra(k) = mean( Dk(triu(true(size(Dk)),1)) , 'omitnan');
    else
        intra(k) = NaN;
    end
    best = inf; bestm = NaN;
    for m = 1:K
        if m==k, continue; end
        Im = find(labels==m);
        if ~isempty(Ik) && ~isempty(Im)
            dm = mean(D(Ik, Im), 'all', 'omitnan');
            if dm < best, best = dm; bestm = m; end
        end
    end
    inter_near(k)  = best;
    nearestCluster(k) = bestm;
end

%% —— margin 单样本检验（星号依据它）——
p_margin = nan(K,1);
margins  = cell(K,1);
for k = 1:K
    Ik = find(labels==k);
    if isempty(Ik) || numel(Ik)<2, p_margin(k)=NaN; margins{k}=nan; continue; end
    % 近本簇
    near_in = arrayfun(@(i) min(D(i, setdiff(Ik,i))), Ik);
    % 近他簇（取使 inter 最小的那簇）
    mstar = nearestCluster(k); Im = find(labels==mstar);
    near_out = min(D(Ik, Im), [], 2);
    mg = near_out - near_in; mg = mg(isfinite(mg));
    margins{k} = mg;
    if isempty(mg)
        p_margin(k)=NaN;
    elseif exist('signrank','file')==2
        p_margin(k) = signrank(mg, 0, 'tail','right');
    else
        Bp=2000; xb = arrayfun(@(~) mean(mg(randsample(numel(mg), numel(mg), true))), 1:Bp);
        p_margin(k) = max(1e-6, mean(xb <= 0));
    end
end

%% 作图（自然顺序 C1..CK）
f = figure('Color','w','Position',[120 120 920 520]); hold on; box on; grid on;

% 颜色
C = [];
if exist('smf_colorlib','file')==2, try, C = smf_colorlib('clusters4_deep'); end, end
if isempty(C), C = lines(max(K,4)); end
if size(C,1) < K, C = repmat(C, ceil(K/size(C,1)), 1); end
C = C(1:K,:);

% 两组柱（同色系：intra=原色；inter=浅色）
X = 1:K; bw=0.36; off=0.18;
b1 = bar(X-off, intra, bw, 'FaceColor','flat','EdgeColor','none');
for i=1:K, b1.CData(i,:) = C(i,:); end
b2 = bar(X+off, inter_near, bw, 'FaceColor','flat','EdgeColor',[0.28 0.28 0.28], 'LineWidth',0.6);
for i=1:K, b2.CData(i,:) = lighten(C(i,:), 0.60); end

% x 轴标签下方加 n（两行标签）
counts = accumarray(max(1,min(K, round(labels(:)))), 1, [K 1])';
xtlbl = arrayfun(@(k) sprintf('C%d\n(n=%d)', k, counts(k)), 1:K, 'UniformOutput', false);
set(gca,'XTick',X, 'XTickLabel', xtlbl, 'TickLabelInterpreter','tex');

ylabel('Mean DTW distance');
legend({'Intra-cluster','Nearest inter-cluster'}, 'Location','northwest','Box','off');
title(sprintf('Fig.2(d)  Intra vs nearest-inter DTW (K=%d)', K), 'FontWeight','bold');
grid(gca,'on'); set(gca,'GridAlpha',0.15);

% —— 在柱顶加：Δ、R 与“更大的星号” —— 
STAR_FS = 13;                 % 星号字号（更醒目）
TXT_FS  = 10;                 % Δ/R 字号
for k = 1:K
    % 效应量
    dlt = inter_near(k) - intra(k);
    rat = inter_near(k) / max(intra(k), eps);
    ytxt = 1.02 * max(intra(k), inter_near(k));
    text(X(k), ytxt, sprintf('\\Delta=%.2f  R=%.2f', dlt, rat), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'Color',[0.20 0.20 0.20], 'FontSize', TXT_FS);

    % 星号（根据 margin 单样本检验）
    stars = p_to_stars(p_margin(k));
    if ~isempty(stars)
        text(X(k), ytxt*1.06, stars, 'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', 'FontWeight','bold', ...
            'Color',[0.15 0.15 0.15], 'FontSize', STAR_FS);
    end
end

% 角标信息
annotation(f,'textbox',[0.60 0.78 0.36 0.16], 'String', ...
    sprintf('K=%d | N=%d | intra=mean upper-tri | inter=min_m mean(D_{k,m}) | stars: margin test', K, numel(labels)), ...
    'EdgeColor',[0.85 0.85 0.85], 'BackgroundColor',[1 1 1], 'FontSize',9);


% 保存
fpng = fullfile(outdir, sprintf('fig2d_main_dtw_K%d.png', K));
fpdf = fullfile(outdir, sprintf('fig2d_main_dtw_K%d.pdf', K));
exportgraphics(f, fpng, 'Resolution', 300);
exportgraphics(f, fpdf);
fprintf('Saved: %s\nSaved: %s\n', fpng, fpdf);

%% 输出
out = struct('labels',labels,'K',K,'intra',intra,'inter_near',inter_near, ...
             'p_margin',p_margin, 'margins',{margins}, 'D',D);
end

% ============================ 辅助函数 ============================
function [labels, K, D] = get_labels_and_D(S)
flds = fieldnames(S);
if numel(flds)==1 && isstruct(S.(flds{1})), S = S.(flds{1}); end
assert(isfield(S,'cluster') && isfield(S.cluster,'labels'), '缺少 S.cluster.labels。');
labels = double(S.cluster.labels(:));
K = max(labels);
if isfield(S.cluster,'K') && ~isempty(S.cluster.K) && S.cluster.K ~= K
    warning('cluster.K=%d 与 labels 最大值 K=%d 不一致；已采用 K=%d（以 labels 为准）。', ...
        S.cluster.K, K, K);
end
assert(isfield(S,'D') && ~isempty(S.D), ['结果文件未包含预计算的 DTW 矩阵 S.D。']);
D = S.D;
N = numel(labels);
assert(all(size(D)==[N N]), 'S.D 尺寸与 labels 不一致：expected %dx%d, got %dx%d.', N, N, size(D,1), size(D,2));
end

function counts = smf_print_cluster_counts(labels, varargin)
p = inputParser;
p.addParameter('Tag','',@(x)ischar(x)||isstring(x));
p.addParameter('ExpectK',[],@(x)isempty(x)||(isscalar(x)&&x==round(x)));
p.parse(varargin{:});
tag = char(p.Results.Tag); Kexp = p.Results.ExpectK;
labels = labels(:);
if any(~isfinite(labels)), warning('labels 中存在 NaN/Inf，将被忽略计数。'); end
L = labels(isfinite(labels));
if isempty(L), fprintf('[%s] labels 为空。\n', tag); counts = []; return; end
K = max(L);
counts = accumarray(max(1,min(K, round(L))), 1, [K 1])';  N = numel(L);
if ~isempty(tag), tinfo = [' [' tag ']']; else, tinfo = ''; end
fprintf('%s Cluster counts (K=%d, N=%d): %s\n', tinfo, K, N, mat2str(counts));
if ~isempty(Kexp) && K ~= Kexp, warning('%s 实际 K=%d 与期望 K=%d 不一致。', tinfo, K, Kexp); end
if sum(counts) ~= N, warning('%s 计数之和(%d) 与 N(%d) 不一致。', tinfo, sum(counts), N); end
for k = 1:K, fprintf('  C%d: %d\n', k, counts(k)); end
end

function c2 = lighten(c, a)
a = max(0, min(1, a));
c2 = (1-a).*c + a.*[1 1 1];
c2 = min(max(c2,0),1);
end

function s = p_to_stars(p)
if isnan(p), s=''; return; end
if p < 1e-3, s='***';
elseif p < 1e-2, s='**';
elseif p < 0.05, s='*';
else, s='';
end
end
