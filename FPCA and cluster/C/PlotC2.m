%% =========== FIG.2c2 — Province × Cluster Heatmap (RUN script, fixed) =============
clear; clc; close all;

% ---------- 配置 ----------
results_mat  = 'smf_fpca_cluster_results.mat';
features_mat = 'features.mat';
outdir       = 'figs_main';
MinPerProvince = 2;
Cmap = parula(256);
if ~exist(outdir,'dir'), mkdir(outdir); end

% ---------- 读取 ----------
S = load(results_mat);
[labels, K] = get_labels_from_results(S);
[variant, province] = get_meta_variant_prov(S, features_mat, numel(labels));

cc = smf_print_cluster_counts(labels, 'Tag','Fig2c2', 'ExpectK',4);
cc_str = strjoin(string(cc), '-');
fname = fullfile(outdir, sprintf('fig2c2_province_heatmap_K%d_[%s]', K, cc_str));



% 省份筛选
prov_all = unique(province,'stable');
keep = false(size(prov_all));
for i=1:numel(prov_all), keep(i) = sum(province==prov_all(i)) >= MinPerProvince; end
prov_use = prov_all(keep);
if isempty(prov_use), warning('无省份满足最小门槛，放宽到全部省份。'); prov_use = prov_all; end
P = numel(prov_use);

% 构建省份×簇计数 O
O = zeros(P, K);
for i=1:P
    idxp = (province==prov_use(i));
    for k=1:K, O(i,k) = sum(idxp & labels==k); end
end
rowsum = sum(O,2); rowsum(rowsum==0)=1;
H = O ./ rowsum;
[Zres, E] = standardized_residuals(O); %#ok<NASGU>

% 行属性
[dom_val, dom_idx] = max(H,[],2);
H_sorted = sort(H,2,'descend');
second_val = zeros(P,1); second_val(:) = H_sorted(:,min(2,K));
ent = row_entropy(H);

% 排序
[~, ord] = sortrows([dom_idx, ent, -second_val, (1:P).'], [1 2 3 4]);
O = O(ord,:); H = H(ord,:); Zres = Zres(ord,:); prov_use = prov_use(ord);
dom_idx = dom_idx(ord); ent = ent(ord);
P = numel(prov_use);

% 右侧变异株细条
VARIANT_ORDER = ["Alpha","Delta","Omicron"];
var_u = VARIANT_ORDER(ismember(VARIANT_ORDER, unique(variant,'stable')));
if isempty(var_u), var_u = unique(variant,'stable'); end
Vbar = zeros(P, numel(var_u));
for i=1:P
    idxp = (province==prov_use(i));
    for j=1:numel(var_u), Vbar(i,j) = sum(idxp & (variant==var_u(j))); end
end
Vbar = Vbar ./ max(1,sum(Vbar,2));

% 颜色
Cclus = smf_colorlib('root');
Cvar  = smf_colorlib('waves');

% ---------- 作图（手动布局，便于拖动/编辑） ----------
f = figure('Color','w','Position',[60 60 1100 760]);

% —— 可调版式参数（单位：normalized）——
margL = 0.07; margR = 0.08; margT = 0.08; margB = 0.08; gapX = 0.02;
wL = 0.14; wM = 0.56; wR = 0.13;           % 左/中/右三块宽度
hAll = 1 - margT - margB;                  % 高度一致

% 计算三块的位置
posL = [margL,                 margB, wL, hAll];
posM = [margL + wL + gapX,     margB, wM, hAll];
posR = [1 - margR - wR,        margB, wR, hAll];

% 创建三块轴（注意：都是可选中/可拖动）
axL = axes('Parent',f,'Units','normalized','Position',posL, ...
    'PickableParts','all','HitTest','on','ActivePositionProperty','position');
axM = axes('Parent',f,'Units','normalized','Position',posM, ...
    'PickableParts','all','HitTest','on','ActivePositionProperty','position');
axR = axes('Parent',f,'Units','normalized','Position',posR, ...
    'PickableParts','all','HitTest','on','ActivePositionProperty','position');

% 左：主导簇徽标
cla(axL); axes_badge(axL, dom_idx, Cclus, P);
yticks(axL, 1:P); yticklabels(axL, cellstr(prov_use));
set(axL,'YDir','normal','XAxisLocation','top');
xlabel(axL,'Dominant'); xlim(axL,[0 1]);
% —— 左侧簇图例（C1–C4）——
hold(axL,'on');
hL = gobjects(1,K);
for k = 1:K
    % 用 NaN 占位的“代表色块”来做图例项，不会影响图面
    hL(k) = patch(axL, nan, nan, Cclus(k,:), ...
                  'EdgeColor','none', 'DisplayName', sprintf('C%d',k), ...
                  'PickableParts','none', 'HitTest','off');
end
legL = legend(axL, hL, arrayfun(@(k)sprintf('C%d',k),1:K,'UniformOutput',false), ...
              'Location','northoutside', 'Orientation','horizontal', ...
              'AutoUpdate','off', 'Box','off');
set(legL, 'Interpreter','none');   % 可在图窗里拖动 legend


% 中：热图
cla(axM); hold(axM,'on');
hImg = imagesc(axM, H, [0 1]); colormap(axM, Cmap);
set(hImg,'PickableParts','all','HitTest','on');     % 热图也可选中
set(axM,'YDir','normal');
xticks(axM,1:K); xticklabels(axM, arrayfun(@(k)sprintf('C%d',k),1:K,'Uni',0));
yticks(axM,1:P); yticklabels(axM,[]); 
cb = colorbar(axM,'Location','eastoutside'); cb.Label.String='Row-normalized share';

% 残差提示点（确保可选）
[rr, cc] = find(abs(Zres)>=2);
for t=1:numel(rr)
    x = cc(t)-0.25; y = rr(t)-0.25;
    plot(axM, x, y, '.', 'Color',[0 0 0], 'MarkerSize',8, ...
        'PickableParts','all','HitTest','on');
end
box(axM,'on'); grid(axM,'off');

% 右：逐省变异株 100% 条
cla(axR); hold(axR,'on'); box(axR,'on');
for i=1:P
    x0 = 0;
    for j=1:size(Vbar,2)
        wj = Vbar(i,j); if wj<=0, continue; end
        ph = patch(axR, [x0 x0+wj x0+wj x0], [i-0.45 i-0.45 i+0.45 i+0.45], ...
              Cvar(1+mod(j-1,size(Cvar,1)),:), 'EdgeColor','none','FaceAlpha',0.95);
        set(ph,'PickableParts','all','HitTest','on');  % 可选可拖
        x0 = x0 + wj;
    end
end
xlim(axR,[0 1]); ylim(axR,[0.5 P+0.5]); set(axR,'YDir','normal');
yticks(axR,1:P); yticklabels(axR,[]); xlabel(axR,'Variants (row-normalized)');
% —— 右侧变异株图例（按 var_u 顺序着色，和右图一致）——
hold(axR,'on');
nVar = numel(var_u);
hR = gobjects(1,nVar);
for j = 1:nVar
    cj = Cvar(1+mod(j-1, size(Cvar,1)), :);
    hR(j) = patch(axR, nan, nan, cj, ...
                  'EdgeColor','none', 'DisplayName', char(var_u(j)), ...
                  'PickableParts','none', 'HitTest','off');
end
legR = legend(axR, hR, cellstr(var_u), ...
              'Location','northoutside', 'Orientation','horizontal', ...
              'AutoUpdate','off', 'Box','off');
set(legR, 'Interpreter','none');   % 可在图窗里拖动 legend


% 同步 y 轴 & 打开编辑
linkaxes([axL,axM,axR],'y');
plotedit on;


% 标题摘要（修正 compose → arrayfun+sprintf）
counts_dom = accumarray(dom_idx,1,[K,1]).';
mean_H = mean(ent);
parts = arrayfun(@(k,c) sprintf('C%d:%d', k, c), 1:K, counts_dom, 'UniformOutput', false);

plotedit on;  % 开启可视化编辑（可选）

ht = annotation(gcf,'textbox',[0.12 0.93 0.76 0.06], ...
    'Units','normalized', ...
    'String', sprintf('Province × Cluster composition — mean entropy = %.2f | Dominant: %s', ...
                      mean_H, strjoin(parts,'  ')), ...
    'EdgeColor','none','HorizontalAlignment','center', ...
    'FontWeight','bold','FontSize',12,'Interpreter','none', ...
    'FitBoxToText','off');


annotation(f,'textbox',[0.12 0.01 0.76 0.04], 'String', ...
    sprintf('Left badges: dominant cluster (C1–C%d).  Right bars: variants %s (row-normalized).  Dots in heatmap: |standardized residual| ≥ 2.', ...
    K, strjoin(cellstr(var_u),'/')), 'EdgeColor','none', 'HorizontalAlignment','center');

% 保存
fname = fullfile(outdir, sprintf('fig2c2_province_heatmap_K%d', K));
print(f,[fname,'.png'],'-dpng','-r300'); savefig(f,[fname,'.fig']);
disp(['Saved: ', fname, '.png & .fig']);

%% ===== 本脚本用到的本地函数 =====
function [labels, K] = get_labels_from_results(S)
    flds = fieldnames(S);
    if numel(flds)==1 && isstruct(S.(flds{1})), S = S.(flds{1}); end
    assert(isfield(S,'cluster') && isfield(S.cluster,'labels'), 'results_mat 缺少 cluster.labels');
    labels = S.cluster.labels(:); K = max(labels);
end
function [variant, province] = get_meta_variant_prov(S, features_mat, N)
    variant = strings(0,1); province = strings(0,1);
    M = []; if isfield(S,'meta'), M=S.meta; elseif isfield(S,'S')&&isfield(S.S,'meta'), M=S.S.meta; end
    pick = @(T, keys) local_pick(T, keys);
    if ~isempty(M)
        variant  = tostr(pick(M,["variant","Variant","strain"]));
        province = tostr(pick(M,["province","Province","region","prov"]));
    end
    if (numel(variant)~=N || numel(province)~=N) && exist(features_mat,'file')
        T = load(features_mat); Tab = pick_widest_table(T);
        if ~isempty(Tab)
            if any(strcmpi(Tab.Properties.VariableNames,'variant')),  variant  = tostr(Tab.variant); end
            if any(strcmpi(Tab.Properties.VariableNames,'province')), province = tostr(Tab.province); end
            if isempty(province) && any(strcmpi(Tab.Properties.VariableNames,'prov')), province = tostr(Tab.prov); end
        end
    end
    if isempty(variant),  variant  = repmat("NA", N,1); end
    if isempty(province), province = repmat("NA", N,1); end
    variant  = tostr(variant); province = tostr(province);
end
function [Z, E] = standardized_residuals(O)
    RT = sum(O,2); CT = sum(O,1); N = sum(RT);
    E = (RT * CT) / max(N, eps);
    Z = (O - E) ./ sqrt(max(E, eps));
end
function H = row_entropy(P)
    P = max(P, eps); P = P ./ sum(P,2);
    H = -sum(P .* log(P), 2);
end
function axes_badge(ax, dom_idx, Cclus, P)
    hold(ax,'on'); xlim(ax,[0 1]); ylim(ax,[0.5 P+0.5]); set(ax,'YDir','normal');
    for i=1:P
        c = Cclus(dom_idx(i),:);
        patch(ax, [0 1 1 0], [i-0.45 i-0.45 i+0.45 i+0.45], c, 'EdgeColor','none');
    end
end
function val = local_pick(T, keys)
    val=[]; if isempty(T), return; end
    if istable(T)
        for i=1:numel(keys), k=char(keys(i));
            if any(strcmpi(T.Properties.VariableNames,k)), val=T.(k); return; end
        end
    elseif isstruct(T)
        f=fieldnames(T);
        for i=1:numel(keys), k=char(keys(i));
            hit=f(strcmpi(f,k)); if ~isempty(hit), val=T.(hit{1}); return; end
        end
    end
end
function Tab = pick_widest_table(T)
    Tab=[]; fn=fieldnames(T); bestW=-1; bestK='';
    for i=1:numel(fn)
        if istable(T.(fn{i})) && width(T.(fn{i}))>bestW
            bestW=width(T.(fn{i})); bestK=fn{i};
        end
    end
    if ~isempty(bestK), Tab=T.(bestK); end
end
function s = tostr(x)
    if isstring(x), s=x; elseif ischar(x), s=string(x);
    elseif iscellstr(x), s=string(x);
    elseif iscell(x), try s=string(x); catch, s=strings(numel(x),1); for i=1:numel(x), if ischar(x{i}), s(i)=string(x{i}); end, end, end
    else, s=string(x); end
end
