function FIG2a_FPCA_BIPLOT()
%% ========================== FIG.2a — 直接运行版 ==========================
% FPCA 形状景观（PC1–PC2）：白底网格 + 原点虚线十字 + 气泡散点（颜色=簇、大小=权重）
% + biplot 箭头（等比例缩短，可单独再缩短最长一根）
%
% 轴装饰（可选，编辑完范围后再执行）：
%   % smf_finalize_ticks(gca, 'HideAxes', true, 'Inside', false, 'OffsetFrac', 0.03, 'Bind', true);
%   % smf_finalize_ticks(gcf, 'HideAxes', true, 'Inside', false, 'OffsetFrac', 0.03, 'Bind', true);
%   % smf_finalize_ticks(gca, 'ShowX', true, 'ShowY', false);
%   % smf_finalize_ticks(gca, 'Remove', true); set(gca,'Visible','on');

clear; clc; close all;

%% 0) 用户配置
results_mat  = 'smf_fpca_cluster_results.mat';
features_mat = 'features.mat';       % 可留空，不存在则跳过箭头
outdir       = 'figs_main';

LABEL_SRC    = 'cluster';            % 'cluster' 或 'byK'
K_TARGET     = 4;                    % 仅当 LABEL_SRC='byK' 时使用

COLOR_MODE   = 'cluster';            % 'cluster' 或 'variant'
PALETTE_KEY  = 'classic4';     % clusters4 / clusters4_alt1/alt2/alt3
SIZE_BY_W    = true;

% —— 点大小增强映射（更显著） ——
SZ_MIN       = 60;                   % 最小点面积
SZ_MAX       = 140;                  % 最大点面积（越大差异越明显）
SIZE_USE_Q   = true;                 % 用分位数拉伸对比（抗极端值）
SIZE_Q_LOW   = 0.05;                 % 下分位
SIZE_Q_HIGH  = 0.95;                 % 上分位
GAMMA_SIZE   = 1.3;                  % γ>1 提升对比（1=线性）

% —— 箭头控制（等比例缩短 & 最长的一根可单独再缩短） ——
ARROW_TOP_N     = 4;                 % 自动选择前 N 个相关性最高的特征
ARROW_NAMES     = {};                % 手动指定：{'delta_early','rho_eps',...}
ARROW_COLOR     = [0.15 0.15 0.15];
ARROW_SCALE     = 1;              % 等比例缩短：0.7=缩到70%，1.0=不缩
LONGEST_SHRINK  = 0.6;              % 仅对"最长那根”再乘此因子（=1 表示不开启）

JITTER_STD   = 0.10;                 % 可加微抖动减少重叠

%% 1) 读取/规整（固定口径：优先用 cluster.labels，K = max(labels)）
S = load(results_mat);
[R, labels, K, src_used] = resolve_labels_fixed(S, LABEL_SRC, K_TARGET);
fprintf('Using labels from: %s | K=max(labels)=%d | N=%d\n', src_used, K, numel(labels));

% 打印四类计数用于核对
smf_print_cluster_counts(labels, 'Tag','Fig2a', 'ExpectK', 4);

w   = R.w;       if isempty(w),   w = ones(numel(labels),1); end
var = R.variant; if isempty(var), var = repmat("NA",numel(labels),1); end
Z   = R.scores(:,1:2);
if JITTER_STD>0, Z = Z + randn(size(Z))*JITTER_STD; end

% 颜色
Cclus = smf_get_colors(PALETTE_KEY, K);  % 统一来自 colorlib
Cvar  = smf_colorlib_safe('variants3');  % 三色（Alpha/Delta/Omicron），循环使用

%% 2) 作图（散点+虚线十字）
if ~exist(outdir,'dir'), mkdir(outdir); end
f = figure('Color','w','Position',[80 80 940 780]); hold on; grid on; box on;

switch lower(COLOR_MODE)
    case 'cluster'
        legNames = arrayfun(@(k)sprintf('C%d',k), 1:K, 'UniformOutput', false);
        for k = 1:K
            idx = (labels==k); if ~any(idx), continue; end
            sz = map_sizes(w(idx), SZ_MIN, SZ_MAX, SIZE_USE_Q, SIZE_Q_LOW, SIZE_Q_HIGH, GAMMA_SIZE, SIZE_BY_W);
            s = scatter(Z(idx,1), Z(idx,2), sz, Cclus(k,:), 'filled', ...
                        'MarkerEdgeColor',[0.15 0.15 0.15], 'LineWidth',0.35);
            s.MarkerFaceAlpha = 0.9; s.MarkerEdgeAlpha = 0.9;
        end
        legend(legNames, 'Location','eastoutside');
    case 'variant'
        uv = unique(var,'stable'); mk = {'o','s','^','d','v','>'};
        for i=1:numel(uv)
            idx = (var==uv(i)); if ~any(idx), continue; end
            ci = 1 + mod(i-1, size(Cvar,1));                 % 循环取色
            sz = map_sizes(w(idx), SZ_MIN, SZ_MAX, SIZE_USE_Q, SIZE_Q_LOW, SIZE_Q_HIGH, GAMMA_SIZE, SIZE_BY_W);
            s = scatter(Z(idx,1), Z(idx,2), sz, Cvar(ci,:), 'filled', ...
                        'Marker', mk{1+mod(i-1,numel(mk))}, ...
                        'MarkerEdgeColor',[0.15 0.15 0.15], 'LineWidth',0.35);
            s.MarkerFaceAlpha = 0.9; s.MarkerEdgeAlpha = 0.9;
        end
        legend(cellstr(uv), 'Location','eastoutside');
    otherwise
        error('COLOR_MODE 仅支持 cluster / variant');
end

% 轴标签（含解释率）
if ~isempty(R.expl_var) && numel(R.expl_var)>=2
    xlabel(sprintf('PC1 (%.1f%%)', 100*R.expl_var(1)));
    ylabel(sprintf('PC2 (%.1f%%)', 100*R.expl_var(2)));
else
    xlabel('PC1'); ylabel('PC2');
end

% 原点虚线十字（先定范围再画）
expand = 0.05;
xlim([min(Z(:,1)) max(Z(:,1))] + [-1 1]*range(Z(:,1))*expand);
ylim([min(Z(:,2)) max(Z(:,2))] + [-1 1]*range(Z(:,2))*expand);
xl = xlim; yl = ylim;
plot([0 0], yl, 'k--', 'LineWidth', 0.9, 'HandleVisibility','off' );
plot(xl, [0 0], 'k--', 'LineWidth', 0.9, 'HandleVisibility','off' );

%% 3) biplot 箭头（来自 features.mat；等比例缩短 + 最长一根再缩短）
if ~isempty(features_mat) && exist(features_mat,'file')
    try
        [Xfeat, names, ids_feat] = extract_features(features_mat);
        % 顺序对齐（若 scores 有 ids）
        if ~isempty(R.ids) && ~isempty(ids_feat)
            [ok, ia, ib] = intersect(ids_feat, R.ids, 'stable');
            if ~isempty(ok)
                Xfeat = Xfeat(ia,:); Zuse = Z(ib,:);
            else
                Zuse = Z;
            end
        else
            Zuse = Z;
        end

        useNames = ARROW_NAMES;
        if isempty(useNames)
            cand = {'delta_early','rho_eps','t_tf2_peak','A_plus','A_minus_eps', ...
                    'kappa_cross','t_cross','t0','tf1_mean_early','tf2_peak'};
            cand = cand(ismember(cand, names)); if isempty(cand), cand = names; end
            useNames = pick_topN_by_corr(Zuse, Xfeat, names, cand, ARROW_TOP_N);
        end
        if ~isempty(useNames)
            draw_biplot_arrows(Zuse, Xfeat, names, useNames, ARROW_COLOR, ARROW_SCALE, LONGEST_SHRINK);
        end
    catch ME
        warning('biplot 箭头阶段跳过：%s', ME.message);
    end
end

title(sprintf('FPCA scatter — source=%s, K=%d', src_used, K));

% —— 保存 PNG + FIG（不再保存 PDF） ——
cc = smf_print_cluster_counts(labels);           % 再打印一次简洁行
cc_str = strjoin(string(cc), '-');               % "25-114-40-31"
fname = fullfile(outdir, sprintf('fig2a_fpca_biplot_%s_K%d_[%s]', lower(src_used), K, cc_str));
if ~exist(outdir,'dir'), mkdir(outdir); end
print(f, [fname,'.png'], '-dpng','-r300');
savefig(f, [fname,'.fig']);
disp(['Saved: ', fname, '.png & .fig']);
end % <<< function 主体结束


%% ============================ 辅助函数们 ============================
function sz = map_sizes(w, szmin, szmax, useQ, ql, qh, gamma_s, enable)
if ~enable, sz = szmin*ones(size(w)); return; end
w = w(:);
if useQ
    lo = quantile(w, ql); hi = quantile(w, qh);
    w = (w - lo) / max(hi - lo, eps);
    w = min(max(w,0),1);
else
    a = min(w); b = max(w); w = (w - a)/max(b-a, eps);
end
if gamma_s ~= 1, w = w.^gamma_s; end
sz = szmin + (szmax - szmin)*w;
end

function [R, labels, K, src_used] = resolve_labels_fixed(S, src, Kreq)
% 显式选择标签来源；若用 cluster，则 K = max(labels)（不信 S.cluster.K）
flds = fieldnames(S);
if numel(flds)==1 && isstruct(S.(flds{1})), S = S.(flds{1}); end
pick = @(T, keys) local_pick(T, keys);

% 常用字段
meta = local_struct_pick(S,'meta'); if isempty(meta), meta = local_struct_pick(S,'metadata'); end
R = struct;
R.scores   = pick(S, ["scores","Z","PC_scores","Scores"]);
R.expl_var = pick(S, ["expl_var","vecr","explained","explained_variance"]);
R.variant  = local_charcell(pick(meta, ["variant","Variant","strain"]));
R.w        = pick(meta, ["weight","w","omega","weights"]);
R.ids      = local_charcell(pick(S, ["ids","ID","wave_id","name","waveID"]));
if ~isempty(R.scores)
    [r,c] = size(R.scores);
    if r < c && r <= 8, R.scores = R.scores.'; end
end

switch lower(src)
    case 'cluster'
        clu = local_struct_pick(S,'cluster');
        assert(~isempty(clu) && isfield(clu,'labels'), 'S.cluster.labels 不存在');
        labels = clu.labels(:);
        K = max(labels);                        % 关键：以实际标签最大值为准
        src_used = 'cluster';
    case 'byk'
        assert(isfield(S,'byK') && ~isempty(S.byK), 'S.byK 不存在');
        KK = arrayfun(@(z) z.K, S.byK);
        pos = find(KK==Kreq,1);
        assert(~isempty(pos), 'byK 中未找到 K=%d', Kreq);
        labels = S.byK(pos).labels(:);
        K = max(labels);                        % 以实际标签最大值为准（防“5类”漏计）
        src_used = sprintf('byK(K=%d)', Kreq);
    otherwise
        error('LABEL_SRC 仅支持 ''cluster'' 或 ''byK''');
end
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
function S2 = local_struct_pick(S,k), S2=[]; if isstruct(S)&&isfield(S,k), S2=S.(k); end, end
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

function [X, names, ids] = extract_features(fmat)
ids = string.empty; X=[]; names=strings(0,1);
T = load(fmat); fn = fieldnames(T); Tab = [];
if numel(fn)==1 && istable(T.(fn{1})), Tab = T.(fn{1});
else
    bestW=-1; bestK='';
    for i=1:numel(fn)
        if istable(T.(fn{i})) && width(T.(fn{i}))>bestW, bestW=width(T.(fn{i}));
            bestK=fn{i};
        end
    end
    if ~isempty(bestK), Tab = T.(bestK); end
end
assert(~isempty(Tab), 'features.mat 未检测到 table 变量。');
idcands = ["wave_id","id","ID","name","waveID"];
hit = intersect(idcands, string(Tab.Properties.VariableNames),'stable');
if ~isempty(hit), ids = string(Tab.(hit(1))); end
isnum = varfun(@(x) isnumeric(x), Tab, 'OutputFormat','uniform');
TabN = Tab(:, isnum); names = string(TabN.Properties.VariableNames); X = table2array(TabN);
end

function chosen = pick_topN_by_corr(Z, X, names, candidates, N)
[tf, loc] = ismember(string(candidates), names); loc = loc(tf);
if isempty(loc), chosen={}; return; end
Xc = X(:, loc);
r1 = corr(Z(:,1), Xc, 'rows','pairwise');
r2 = corr(Z(:,2), Xc, 'rows','pairwise');
mag = sqrt(r1.^2 + r2.^2);
[~, ord] = sort(mag, 'descend');
ord = ord(1:min(N, numel(ord)));
chosen = cellstr(names(loc(ord)));
end

function draw_biplot_arrows(Z, X, names, useNames, arrowColor, arrowScale, longestShrink)
[tf, loc] = ismember(string(useNames), names); loc = loc(tf);
if isempty(loc), return; end
Xc = X(:, loc);
r1 = corr(Z(:,1), Xc, 'rows','pairwise');
r2 = corr(Z(:,2), Xc, 'rows','pairwise');
V  = [r1(:), r2(:)];

% 归一到画幅，再整体缩放
xl = xlim; yl = ylim; base = 0.80 * min(range(xl), range(yl));
nrm = sqrt(sum(V.^2,2)); if any(nrm>0), V = V ./ max(nrm); end
V = V * base * arrowScale;

% 仅对最长的一根再缩短（可选）
if longestShrink ~= 1
    lens = vecnorm(V,2,2);
    [~, imax] = max(lens);
    V(imax,:) = V(imax,:) * longestShrink;
end

for i=1:size(V,1)
    quiver(0,0, V(i,1), V(i,2), 0, 'Color', arrowColor, 'LineWidth',1.2, ...
           'MaxHeadSize',0.25, 'HandleVisibility','off' );
    tx = V(i,1) * 1.05; ty = V(i,2) * 1.05;
    text(tx, ty, useNames{i}, 'FontSize', 10, 'Color', arrowColor, ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'BackgroundColor',[1 1 1],'Margin',1);
end
end

function counts = smf_print_cluster_counts(labels, varargin)
% 打印聚类计数，便于快速核对是否一致
p = inputParser;
p.addParameter('Tag','',@(x)ischar(x)||isstring(x));
p.addParameter('ExpectK',[],@(x)isempty(x)||(isscalar(x)&&x==round(x)));
p.parse(varargin{:});
tag = char(p.Results.Tag);
Kexp = p.Results.ExpectK;

labels = labels(:);
if any(~isfinite(labels))
    warning('labels 中存在 NaN/Inf，将被忽略计数。');
end
L = labels(isfinite(labels));
if isempty(L)
    fprintf('[%s] labels 为空。\n', tag); counts = []; return;
end
K = max(L);
counts = accumarray(max(1,min(K, round(L))), 1, [K 1])';
N = numel(L);

if ~isempty(tag), tinfo = ['[' tag '] ']; else, tinfo = ''; end
fprintf('%sCluster counts (K=%d, N=%d): %s\n', tinfo, K, N, mat2str(counts));
for k = 1:K, fprintf('  C%d: %d\n', k, counts(k)); end
if ~isempty(Kexp) && K ~= Kexp
    warning('%s实际 K=%d 与期望 K=%d 不一致。', tinfo, K, Kexp);
end
if sum(counts) ~= N
    warning('%s计数之和(%d) 与 N(%d) 不一致。', tinfo, sum(counts), N);
end
end

function C = smf_get_colors(palette_key, K)
% 统一从 smf_colorlib 取色；若库缺失或不足 K 行则循环补齐
Cbase = smf_colorlib_safe(palette_key);
if size(Cbase,1) >= K
    C = Cbase(1:K,:);
else
    rep = ceil(K/size(Cbase,1));
    C = repmat(Cbase, rep, 1);
    C = C(1:K,:);
end
end

function C = smf_colorlib_safe(which)
% 安全调用：若用户未提供 smf_colorlib.m，则给出默认配色
if exist('smf_colorlib','file')==2
    C = smf_colorlib(which);
    return;
end
% 默认：clusters4_alt2 / variants3 的兜底
switch lower(which)
    case {'clusters4','clusters4_alt1','clusters4_alt2','clusters4_alt3'}
        C = [0.10 0.28 0.46; 0.78 0.64 0.00; 0.13 0.55 0.28; 0.75 0.11 0.20]; % alt2
    case 'variants3'
        C = [0.12 0.47 0.71; 1.00 0.50 0.05; 0.17 0.63 0.17];
    otherwise
        C = lines(4);
end
end
