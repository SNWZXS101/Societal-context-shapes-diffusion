function FIG2e_Representatives_FromMain()
%% Fig.2(e) — Representatives per cluster + prototype（精简独立版）
% 读取 smf_fpca_cluster_results.mat，按你的主脚本(e)口径：
%   - 每簇：medoid（若有） + 按权重 Top 的成员（总数≤REPS_PER_CLUSTER）
%   - 与聚类原型叠加显示
% 颜色与计数：smf_colorlib / smf_print_cluster_counts（与你 Fig2d 一致）

clear; clc; close all;

%% 配置
results_mat      = 'smf_fpca_cluster_results.mat';
features_mat     = 'features.mat';
outdir           = 'figs_main';
REPS_PER_CLUSTER = 3;                % 每簇展示的成员条数（包含 medoid）
PALETTE_KEY      = 'clusters4_deep';

%% 读取结果并统一字段到 R（复用你主脚本的口径）
S = load(results_mat);
R = resolve_result_struct_local(S);          % 保证 R.u, R.F_norm, R.scores 等
[labels, proto, K] = choose_labels_prototypes_local(S, R);  % K=实际 labels 最大值
if ~exist(outdir,'dir'), mkdir(outdir); end

% —— 原型朝向规范（与主脚本一致：proto: M×K 且 M==numel(R.u)）——
if ~isempty(proto) && ~isempty(R.u)
    [rp,cp] = size(proto);
    if rp == numel(R.u)
        % ok
    elseif cp == numel(R.u)
        proto = proto.';  % to M x K
    end
end

% —— 权重（优先 results.meta.weight；否则 features 的 w_*；长度不符则全 1）——
w = meta_from_results_or_features_local(S, features_mat);
N = size(R.F_norm,2);
if ~isvector(w) || numel(w)~=N
    w = ones(N,1);
else
    w = double(w(:));
    bad = ~isfinite(w); if any(bad), w(bad) = median(w(~bad)); end
    w(w<0)=0;
end

% 打印计数（你的函数）
if exist('smf_print_cluster_counts','file')==2
    smf_print_cluster_counts(labels, 'Tag','Fig2e', 'ExpectK', K);
end

% 颜色（你的 colorlib；不足则循环补齐）
if exist('smf_colorlib','file')==2
    try
        C = smf_colorlib(PALETTE_KEY);
    catch
        C = lines(max(4,K));
    end
else
    C = lines(max(4,K));
end
if size(C,1) < K
    C = repmat(C, ceil(K/size(C,1)), 1);
end
C = C(1:K,:);

%% 代表索引（完全沿用你主脚本逻辑）
Reps = get_representatives_local(S, labels, REPS_PER_CLUSTER, w);

%% 作图
f = figure('Color','w','Position',[80 40 1200 650]);
tiledlayout(ceil(K/2),2,'TileSpacing','compact','Padding','compact');

for k = 1:K
    nexttile; hold on; grid on; box on;

    % 成员
    reps = Reps{k};
    for r = 1:numel(reps)
        j = reps(r);
        if j>=1 && j<=size(R.F_norm,2)
            plot(R.u, R.F_norm(:,j), '-', 'Color',[C(k,:) 0.85], 'LineWidth',1.4);
        end
    end

    % 原型（若无则不画；与主脚本保持一致）
    if ~isempty(proto) && size(proto,2)>=k
        plot(R.u, proto(:,k), '-', 'Color',C(k,:), 'LineWidth',2.6);
    end

    % 小装饰
    idx = find(labels==k);
    title(sprintf('C%d  |  n=%d  |  reps=%d', k, numel(idx), numel(reps)));
    xlim([min(R.u) max(R.u)]);
    xlabel('Normalized time (0–1)'); ylabel('Normalized amplitude');
end

sgtitle(sprintf('Fig.2(e) Representatives vs prototype (K=%d)', K), 'FontWeight','bold');

fname = fullfile(outdir, sprintf('fig2e_representatives_K%d',K));
print(f,[fname,'.png'],'-dpng','-r300'); savefig(f,[fname,'.fig']);
disp('Fig.2e done.');
end

%% ============================ 辅助函数们（与主脚本口径一致） ============================
function R = resolve_result_struct_local(S)
R = struct;
flds = fieldnames(S);
if numel(flds)==1 && isstruct(S.(flds{1})), S = S.(flds{1}); end

pick = @(T, keys) local_pick(T, keys);

R.u        = pick(S, ["u","grid","time_grid","U"]);
R.F_norm   = pick(S, ["F_norm","F","curves","X","Y"]);
R.mean_f   = pick(S, ["mean_f","mu","mean_curve","MeanF"]);
R.Phi      = pick(S, ["Phi","loadings","eigenfunctions","PC_funcs"]);
R.scores   = pick(S, ["scores","Z","PC_scores","Scores"]);
R.expl_var = pick(S, ["expl_var","vecr","explained","explained_variance"]);
R.proto    = pick(struct_pick(S,"cluster"), ["prototypes","proto","centers","DBA"]);
R.labels   = pick(struct_pick(S,"cluster"), ["labels","idx","cluster_id","assignments"]);
M = struct_pick(S, "meta"); if isempty(M), M = struct_pick(S, "metadata"); end
R.province = charcell(pick(M, ["province","prov","Province"]));
R.variant  = charcell(pick(M, ["variant","Variant","strain"]));
R.w        = pick(M, ["weight","w","omega","weights"]);

% 必要校验与朝向
assert(~isempty(R.u) && ~isempty(R.F_norm), 'Missing u/F_norm in results.');
if ~isempty(R.scores)
    [r,c] = size(R.scores);
    if r < c && r <= 8, R.scores = R.scores.'; end
end
if ~isempty(R.u) && ~isempty(R.F_norm)
    [r,c] = size(R.F_norm);
    if r ~= numel(R.u) && c == numel(R.u)
        R.F_norm = R.F_norm.'; % M x N，与主脚本一致
    end
end
if ~isempty(R.mean_f)
    if isrow(R.mean_f), R.mean_f = R.mean_f.'; end
end
end

function [labels, proto, K] = choose_labels_prototypes_local(S, R)
% 与主脚本：优先 cluster.labels，K=labels 最大值；proto 用 S.cluster.prototypes
labels = []; proto = []; K = [];
flds = fieldnames(S);
if numel(flds)==1 && isstruct(S.(flds{1})), S = S.(flds{1}); end
if isfield(S,'cluster') && isstruct(S.cluster)
    labels = S.cluster.labels;
    proto  = S.cluster.prototypes;
    K = max(labels);
else
    error('No cluster info in results.');
end
% 原型朝向在主流程中再处理（依 R.u 转置）
end

function Reps = get_representatives_local(S, labels, kPerCluster, w)
K = max(labels);
Reps = cell(K,1);
haveReps = (isfield(S,'cluster') && isfield(S.cluster,'representatives') && ~isempty(S.cluster.representatives));
haveMed  = (isfield(S,'cluster') && isfield(S.cluster,'medoids') && ~isempty(S.cluster.medoids));

if haveReps
    Rcell = S.cluster.representatives;
    if iscell(Rcell)
        for k=1:K
            idx = Rcell{k}; if isempty(idx), idx = []; end
            Reps{k} = idx(:).';
        end
    elseif isnumeric(Rcell)
        for k=1:K, Reps{k} = Rcell(:,k).'; end
    else
        for k=1:K, Reps{k} = []; end
    end
elseif haveMed
    % 用 medoid + 权重 Top（排除已选）
    for k=1:K
        idx = find(labels==k).';
        cand = [];
        m = S.cluster.medoids(min(k,numel(S.cluster.medoids)));
        if isstruct(m), m = m.medoids; end
        m = round(m);
        if isfinite(m) && m>=1 && m<=numel(labels), cand = [cand, m]; end %#ok<AGROW>
        idxRemain = setdiff(idx, cand, 'stable');
        if ~isempty(idxRemain)
            [~, ord] = sort(w(idxRemain),'descend');
            take = min(kPerCluster - numel(cand), numel(ord));
            if take>0, cand = [cand, idxRemain(ord(1:take))]; end %#ok<AGROW>
        end
        Reps{k} = cand;
    end
else
    % 仅按权重 Top
    for k=1:K
        idx = find(labels==k);
        if isempty(idx), Reps{k} = []; continue; end
        [~,ord] = sort(w(idx),'descend');
        Reps{k} = idx(ord(1:min(kPerCluster,numel(ord))));
    end
end
end

function w = meta_from_results_or_features_local(S, features_mat)
% 返回单输出 w；优先 results.mat 的 meta.weight，其次 features.mat 的 w_tail_14/w_tail/weight
w = [];
% results 中 meta
M = [];
if isstruct(S)
    if isfield(S,'meta')
        M = S.meta;
    elseif isfield(S,'S') && isfield(S.S,'meta')
        M = S.S.meta;
    end
end
if ~isempty(M)
    tmp = local_pick(M, ["weight","w","omega","weights"]);
    if ~isempty(tmp), w = tmp(:); end
end
% features 兜底
if (isempty(w) || ~isvector(w)) && exist(features_mat,'file')
    T = load(features_mat); fn = fieldnames(T); Tab=[];
    if isfield(T,'Features') && istable(T.Features), Tab=T.Features; end
    if isempty(Tab)
        bestW=-1; bestK=''; 
        for i=1:numel(fn)
            ti = T.(fn{i});
            if istable(ti) && width(ti)>bestW, bestW=width(ti); bestK=fn{i}; end
        end
        if ~isempty(bestK), Tab = T.(bestK); end
    end
    if ~isempty(Tab)
        vnames = string(Tab.Properties.VariableNames);
        key = intersect(lower(vnames), {'w_tail_14','w_tail','weight'}, 'stable');
        if ~isempty(key)
            % 找到原名（避免大小写）
            kreal = vnames{strcmpi(vnames, key{1})};
            w = Tab.(kreal); w = w(:);
        end
    end
end
if ~isempty(w), w = double(w(:)); end
end

%% ====== 通用小工具（与你主脚本相同口径） ======
function val = local_pick(T, names)
val = [];
if isempty(T), return; end
if istable(T)
    vn = T.Properties.VariableNames;
    for i=1:numel(names)
        k = char(names(i));
        if any(strcmpi(vn, k)), val = T.(k); return; end
    end
elseif isstruct(T)
    f = fieldnames(T);
    for i=1:numel(names)
        k = char(names(i));
        hit = f(strcmpi(f,k));
        if ~isempty(hit), val = T.(hit{1}); return; end
    end
end
end

function S2 = struct_pick(S,k)
S2 = [];
if isstruct(S) && isfield(S,k), S2 = S.(k); end
end

function out = charcell(x)
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
