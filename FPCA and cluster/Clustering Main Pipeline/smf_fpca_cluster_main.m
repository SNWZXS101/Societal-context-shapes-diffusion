function smf_fpca_cluster_main(varargin)
% SMF_FPCA_CLUSTER_MAIN  (features-only, full version)
% End-to-end pipeline for SMF shape comparison with:
% - load ONLY from features.mat (authority source)
% - weights = features.w_tail_14
% - A+ normalization on common grid (strict check)
% - weighted FPCA
% - DTW + K-medoids + DBA (scan K)
% - save per-K results (byK), distance matrix D, kNN
% - "silhouette_unique" (de-duplicated within-cluster by DTW-eps)
% - stability across normalizations: f~1 / f~2 (if present in features)
%
% Name-Value parameters (key ones):
%   'features_mat'  (char)   : path to features.mat   [default 'features.mat']
%   'M'             (int)    : grid size              [200]
%   'K_grid'        (vector) : K list to scan         [3:8]
%   'K_default'     (int)    : default K to select    [5]
%   'DTW_Window'    (scalar) : Sakoe–Chiba band       [0.2]
%   'UseParallel'   (bool)   : parallel if available  [false]
%   'Nmax'          (int)    : take first N rows      []
%   'indices'       (vector) : take specific rows     []
%   'Unique_eps'    (scalar) : eps for silhouette_unique threshold [1e-6]
%
% Optional hints for alt channels (if your features table uses special col names):
%   'AltF1Names'    (cellstr): candidate column names for f~1 series
%   'AltF2Names'    (cellstr): candidate column names for f~2 series
%   'TimeNames'     (cellstr): candidate column names for time series

% ------------- params -------------
p = inputParser;
addParameter(p,'features_mat','features.mat',@ischar);
addParameter(p,'M',200,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'K_grid',3:8,@(x)isnumeric(x)&&isvector(x));
addParameter(p,'K_default',3,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'dba_iter',15,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'knn_k',6,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'DTW_Window',0.1,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'UseParallel',false,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'Nmax',[],@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'indices',[],@(x)isnumeric(x)&&isvector(x));
addParameter(p,'Unique_eps',1e-6,@(x)isnumeric(x)&&isscalar(x));

addParameter(p,'MultipartAsSeparate', true, @(x)islogical(x)||ismember(x,[0 1]));   % 是否将多波单独成类
addParameter(p,'MultipartLabelAs','K+1', @(x) any(strcmpi(x,{'K+1','0'})));         % 多波标签：K+1 或 0
addParameter(p,'NormMode','Aplus', @(x) any(strcmpi(x,{'Aplus','Fmax'})));


% optional alt-channel name hints
addParameter(p,'AltF1Names',{'f_tilde1','ftilde1','f1','f_beta','f_betaN'},@(x)iscellstr(x)||isstring(x));
addParameter(p,'AltF2Names',{'f_tilde2','ftilde2','f2','f_betaSI','f_over_betaSI'},@(x)iscellstr(x)||isstring(x));
addParameter(p,'TimeNames', {'T','t','time','x','u','grid'},@(x)iscellstr(x)||isstring(x));

parse(p,varargin{:});
args = p.Results;


% ------------- load (features only) -------------
[W, Tfeat] = load_from_features(args.features_mat);

N = numel(W.f);

% indices subset first (preferred)
if ~isempty(args.indices)
    I = unique(args.indices(:));
    I(I<1 | I>N) = [];
    assert(~isempty(I), 'indices empty after sanitization.');
    W.f = W.f(I); W.T = W.T(I);
    W.province = W.province(I); W.variant = W.variant(I);
    W.wave_id  = W.wave_id(I);  W.wave_idx = W.wave_idx(I);
    Tfeat = Tfeat(I,:);
    N = numel(W.f);
end
% Nmax subset (secondary)
if ~isempty(args.Nmax) && args.Nmax < N
    I = 1:args.Nmax;
    W.f = W.f(I); W.T = W.T(I);
    W.province = W.province(I); W.variant = W.variant(I);
    W.wave_id  = W.wave_id(I);  W.wave_idx = W.wave_idx(I);
    Tfeat = Tfeat(I,:); N = numel(W.f);
end
fprintf('Loaded %d waves from %s (features only)\n', N, args.features_mat);

% === 拆分单波 / 多波 ===
if ~ismember('flag_multipart', Tfeat.Properties.VariableNames)
    error('Features 缺少 flag_multipart 列（请确保 Excel 覆盖已生效）');
end
idxM = (Tfeat.flag_multipart == 1);   % 多波
idxS = ~idxM;                         % 单波
fprintf('Split: single=%d, multipart=%d (total=%d)\n', sum(idxS), sum(idxM), N);


% ------------- weights = w_tail_14 -------------
[omega, used_col] = derive_weights(Tfeat);
assert(numel(omega)==N, 'w_tail_14 length mismatch with number of waves.');
omega = double(omega(:));

% ------------- grid + normalization (Aplus | Fmax) -------------
u = linspace(0,1,args.M)';
switch lower(args.NormMode)
    case 'aplus'
        [F_Ap, scale_after] = resample_Aplus_normalize(W.f, W.T, u, W.wave_id);
    case 'fmax'
        [F_Ap, scale_after] = resample_fmax_normalize(W.f, W.T, u, W.wave_id);
    otherwise
        error('Unknown NormMode: %s', args.NormMode);
end
fprintf('Scale (%s) after (grid) => min=%.4f, max=%.4f, mean=%.4f\n', ...
        args.NormMode, min(scale_after), max(scale_after), mean(scale_after));


% ------------- FPCA -------------
if args.MultipartAsSeparate
    % 仅用单波的曲线与权重来拟合主成分
    [Phi, Zs, expl_var, mu] = weighted_fpca(F_Ap(:, idxS), omega(idxS));
    % 将全体样本投影到同一基（便于统一可视化/保存）
    scores = ((F_Ap - mu)') * Phi;   % N x nPC
else
    [Phi, scores, expl_var, mu] = weighted_fpca(F_Ap, omega);
end

% ------------- DTW distances -------------
fprintf('Computing pairwise DTW distances (N=%d, M=%d)...\n', N, args.M);
D = pairwise_dtw(F_Ap, 'Window', args.DTW_Window, 'UseParallel', args.UseParallel);

% ------------- scan K -------------
Klist = args.K_grid(:)';
qual = struct('K_grid',Klist, ...
              'silhouette',nan(size(Klist)), ...
              'silhouette_unique',nan(size(Klist)), ...
              'gap',nan(size(Klist)));
labels_by_K = cell(size(Klist));
protos_by_K = cell(size(Klist));
medoids_by_K= cell(size(Klist));

if ~args.MultipartAsSeparate
    % ====== 原逻辑：全体一起聚类 ======
    for ii = 1:numel(Klist)
        K = Klist(ii);
        [idx, medoid_idx] = kmedoids_precomputed(D, K);
        protos = dba_prototypes(F_Ap, idx, K, omega, args.dba_iter);
        silw = weighted_silhouette(D, idx, omega);
        qual.silhouette(ii) = mean(silw(~isnan(silw)));
        qual.silhouette_unique(ii) = silhouette_unique_fn(D, idx, args.Unique_eps);
        qual.gap(ii) = gap_like(D, idx, 50);
        labels_by_K{ii} = idx;  protos_by_K{ii} = protos;  medoids_by_K{ii}= medoid_idx;
        fprintf('K=%d | silhouette=%.3f | silhouette_unique=%.3f\n', ...
                K, qual.silhouette(ii), qual.silhouette_unique(ii));
    end
else
    % ====== 新逻辑：只对“单波子集”聚类，多波单独成类 ======
    DS = D(idxS, idxS);                  % 单波子集距离
    omegaS = omega(idxS);
    S_inds = find(idxS);                 % 单波的全局索引（映射 medoid 用）

    for ii = 1:numel(Klist)
        K = Klist(ii);

        % 在单波子集上 K-medoids
        [idxS_labels, medoid_idx_S] = kmedoids_precomputed(DS, K);

        % 质量评估（只在单波子集上）
        silwS = weighted_silhouette(DS, idxS_labels, omegaS);
        qual.silhouette(ii) = mean(silwS(~isnan(silwS)));
        qual.silhouette_unique(ii) = silhouette_unique_fn(DS, idxS_labels, args.Unique_eps);
        qual.gap(ii) = gap_like(DS, idxS_labels, 50);

        % 构造“全体”的标签（多波并回）
        lab_full = zeros(N,1);
        lab_full(idxS) = idxS_labels;
        if strcmpi(args.MultipartLabelAs,'K+1')
            lab_full(idxM) = K + 1;
        else
            lab_full(idxM) = 0;   % 0 表示未聚类/多波
        end

        % 原型：单波各簇用 DBA，多波（若为K+1）可选做一个 DBA 原型
        protosS = dba_prototypes(F_Ap(:, idxS), idxS_labels, K, omegaS, args.dba_iter);   % M x K
        if any(idxM) && strcmpi(args.MultipartLabelAs,'K+1')
            protoM = dba_barycenter(F_Ap(:, idxM), omega(idxM), args.dba_iter);           % M x 1
            protos_full = [protosS, protoM];
        else
            protos_full = protosS;
        end

        % medoid：把子集 medoid 映射回全局索引；多波给一个代表（可选）
        medoids_global = S_inds(medoid_idx_S(:));
        if any(idxM) && strcmpi(args.MultipartLabelAs,'K+1')
            medoids_global = [medoids_global; find(idxM,1,'first')];
        end

        labels_by_K{ii}  = lab_full;
        protos_by_K{ii}  = protos_full;
        medoids_by_K{ii} = medoids_global;

        fprintf('K=%d | silhouette=%.3f | silhouette_unique=%.3f (single-wave)\n', ...
                K, qual.silhouette(ii), qual.silhouette_unique(ii));
    end
end


% ------------- choose default K -------------
Kstar = args.K_default;
kpos = find(Klist==Kstar,1);
if isempty(kpos), [~,kpos] = max(qual.silhouette); Kstar = Klist(kpos); end
labels     = labels_by_K{kpos};
prototypes = protos_by_K{kpos};
medoids    = medoids_by_K{kpos};
repr       = pick_representatives(F_Ap, labels, prototypes, 6);

% ------------- kNN for advanced viz -------------
knn_idx = build_knn_from_D(D, args.knn_k);

% ------------- stability across normalizations (ftilde1/ftilde2 if present) -------------
stab = struct('ARI',[],'VI',[],'labels',{{}});
altnorm = struct();

try
    [F_t1, has_f1] = resample_alt_from_features(Tfeat, args.TimeNames, args.AltF1Names, u);
    [F_t2, has_f2] = resample_alt_from_features(Tfeat, args.TimeNames, args.AltF2Names, u);

    label_sets = {labels};
    names = {"Aplus"};

    if has_f1
        D1 = pairwise_dtw(F_t1, 'Window', args.DTW_Window, 'UseParallel', false);
        [idx1, ~] = kmedoids_precomputed(D1, Kstar);
        [~, Z1] = weighted_fpca(F_t1, omega); 
        altnorm.ftilde1 = struct('scores',Z1(:,1:2),'cluster',struct('labels',idx1));
        label_sets{end+1} = idx1; names{end+1} = 'ftilde1';
    end
    if has_f2
        D2 = pairwise_dtw(F_t2, 'Window', args.DTW_Window, 'UseParallel', false);
        [idx2, ~] = kmedoids_precomputed(D2, Kstar);
        [~, Z2] = weighted_fpca(F_t2, omega); %#ok<ASGLU>
        altnorm.ftilde2 = struct('scores',Z2(:,1:2),'cluster',struct('labels',idx2));
        label_sets{end+1} = idx2; names{end+1} = 'ftilde2';
    end

    m = numel(label_sets);
    ARI = nan(m); VI = nan(m);
    for i=1:m
        for j=1:m
            ARI(i,j) = ari(label_sets{i}, label_sets{j});
            VI(i,j)  = vi(label_sets{i},  label_sets{j});
        end
    end
    stab.ARI = ARI; stab.VI = VI; stab.labels = names;

catch ME
    warning('SMF:StabilitySkipped', ...
        'Stability across normalizations skipped (%s): %s', ...
        ME.identifier, ME.message);
end

% ------------- pack & save (with byK, D) -------------
cluster = struct('K',Kstar,'labels',labels,'prototypes',prototypes, ...
                 'representatives',{repr}, 'medoids',medoids);
quality = struct('K_grid',qual.K_grid,'silhouette',qual.silhouette, ...
                 'silhouette_unique',qual.silhouette_unique,'gap',qual.gap);

byK = struct('K',[],'labels',[],'prototypes',[],'medoids',[], ...
             'silhouette',[],'silhouette_unique',[],'gap',[]);
byK = repmat(byK, numel(Klist), 1);
for ii = 1:numel(Klist)
    byK(ii).K                 = Klist(ii);
    byK(ii).labels            = labels_by_K{ii};
    byK(ii).prototypes        = protos_by_K{ii};
    byK(ii).medoids           = medoids_by_K{ii};
    byK(ii).silhouette        = qual.silhouette(ii);
    byK(ii).silhouette_unique = qual.silhouette_unique(ii);
    byK(ii).gap               = qual.gap(ii);
end

meta = struct();
meta.province = W.province;
meta.variant  = W.variant;
meta.wave_id  = W.wave_id;
meta.wave_idx = W.wave_idx;
meta.weight   = omega(:);
meta.weight_source = used_col;

norminfo = struct('mode',args.NormMode,'scale',scale_after);


Results = struct();
Results.u          = u;
Results.F_norm     = F_Ap;
Results.mean_f     = mu;
Results.Phi        = Phi;
Results.scores     = scores;
Results.expl_var   = expl_var;
Results.cluster    = cluster;
Results.byK        = byK;
Results.cluster_quality = quality;
Results.D          = D;
Results.knn        = struct('idx',knn_idx);
Results.stability  = stab;
Results.altnorm    = altnorm;
Results.meta       = meta;
if strcmpi(args.NormMode,'Aplus')
    Results.Aplus_after = scale_after;   % 兼容旧字段名
else
    Results.Aplus_after = [];            % 或干脆不写这个字段
end
Results.norm = norminfo;                 % 新字段：规范记录当前归一方式




save('smf_fpca_cluster_results.mat','-struct','Results','-v7.3');
fprintf('Saved: smf_fpca_cluster_results.mat\n');
end

% ================== local helpers (features-only) ==================

function s_mean = silhouette_unique_fn(Dfull, labels, eps0)
% 去重后 silhouette（只在簇内合并 DTW<=eps0 的近重复）
N = size(Dfull,1); %#ok<NASGU>
K = max(labels);
rep  = [];  Lg = [];  wg = [];
for k = 1:K
    idxk = find(labels==k);
    if isempty(idxk), continue; end
    Dk = Dfull(idxk, idxk);
    A = (Dk <= eps0); A = A|A.'; A(1:size(A,1)+1:end) = true;
    G = graph(A,'upper'); bins = conncomp(G);
    for g = 1:max(bins)
        members = idxk(bins==g);
        rep(end+1,1) = members(1);
        Lg(end+1,1)  = k;
        wg(end+1,1)  = numel(members);
    end
end
Dg = Dfull(rep, rep);
s  = weighted_silhouette(Dg, Lg, wg);
s(~isfinite(s)) = 0;
s_mean = mean(s);
end

function [F_alt, has] = resample_alt_from_features(Tfeat, timeNameCands, altNameCands, u)
% 从 features 表中查找“备用通道”的曲线列（cell of numeric），插值到 u
has = false; F_alt = [];
vn = string(Tfeat.Properties.VariableNames);
% time col
tcol = "";
for k = 1:numel(timeNameCands)
    hit = find(strcmpi(vn, timeNameCands{k}), 1);
    if ~isempty(hit) && iscell(Tfeat.(hit)) && ~isempty(Tfeat.(hit)) && isnumeric(Tfeat.(hit){1})
        tcol = vn(hit); break;
    end
end
% alt f col
fcol = "";
for k = 1:numel(altNameCands)
    hit = find(strcmpi(vn, altNameCands{k}), 1);
    if ~isempty(hit) && iscell(Tfeat.(hit)) && ~isempty(Tfeat.(hit)) && isnumeric(Tfeat.(hit){1})
        fcol = vn(hit); break;
    end
end
if strlength(fcol)==0 || strlength(tcol)==0
    return;
end
N = height(Tfeat); M = numel(u);
F_alt = nan(M,N);
for i=1:N
    f = double(Tfeat.(fcol){i}(:));
    t = double(Tfeat.(tcol){i}(:));
    if numel(t)>=2 && (max(t)>min(t))
        t0 = (t - min(t)) / (max(t) - min(t));
        F_alt(:,i) = interp1(t0, f, u, 'pchip','extrap');
    else
        F_alt(:,i) = repmat(mean(f,'omitnan'), M, 1);
    end
end
has = true;
end
