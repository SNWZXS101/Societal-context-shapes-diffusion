%% PLOT_EXTENDED_FIGURES — 自包含升级版
% (a) Scree plot (explained variance)
% (b) Cluster quality across K (silhouette / silhouette_unique / gap-like)
% (c) Cluster temporal prevalence (若 meta.date 可解析)

clear; clc; close all;

%% ===== 路径与参数 =====
results_mat = 'smf_fpca_cluster_results.mat';
outdir      = 'figs_ext';
ensure_dir(outdir);

Win         = 0.20;   % 若需重算 DTW
Gap_B       = 200;    % GAP 统计自举次数
rng(1);               % 固定种子以保证 GAP 可重复

%% ===== 读取与规整 =====
S = load(results_mat);
R = resolve_result_struct_ext(S);    % 自带 expl_var & cluster_quality

%% >>> 若 cluster_quality 缺失或不完整，自动补算 <<<
needCalc = isempty(R.K_grid) || isempty(R.silhouette);
if needCalc
    fprintf('[Info] cluster_quality 缺失，自动补算…\n');
    [Klist, sil, silU, gap] = calc_quality_byK(S, Win, Gap_B);
    R.K_grid            = Klist;
    R.silhouette        = sil;
    R.silhouette_unique = silU;
    R.gap               = gap;
else
    % 转成列向量，方便后续画图
    R.silhouette        = R.silhouette(:);
    R.silhouette_unique = R.silhouette_unique(:);
    R.gap               = R.gap(:);
end

%% (a) Scree — 解释率
if ~isempty(R.expl_var)
    figure('Color','w','Position',[100 100 700 480]);
    bar(100*R.expl_var(:), 'FaceColor',[0.20 0.50 0.80]); hold on;
    plot(cumsum(100*R.expl_var(:)),'-k','LineWidth',1.5);
    grid on; xlabel('PC'); ylabel('Explained variance (%)');
    title('FPCA explained variance and cumulative');
    print(fullfile(outdir,'ext_scree.png'),'-dpng','-r300');
    print(fullfile(outdir,'ext_scree.pdf'),'-dpdf','-r300');
end

%% (b) Cluster quality vs K
if ~isempty(R.K_grid)
    figure('Color','w','Position',[100 100 900 480]);
    tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    % (左) silhouette
    nexttile; hold on;
    if ~isempty(R.silhouette)
        plot(R.K_grid, R.silhouette,'-o','LineWidth',1.8);
    end
    if ~isempty(R.silhouette_unique)
        plot(R.K_grid, R.silhouette_unique,'-s','LineWidth',1.8);
    end
    grid on; xlabel('K'); ylabel('Score'); title('Silhouette vs K');
    if ~isempty(R.silhouette_unique)
        legend({'silhouette','silhouette unique'},'Location','best');
    end

    % (右) gap
    nexttile;
    if ~isempty(R.gap)
        plot(R.K_grid, R.gap,'-d','LineWidth',1.8);
        grid on; xlabel('K'); ylabel('Gap-like'); title('Gap vs K');
    else
        axis off; text(0.1,0.5,'No gap statistics','FontSize',15);
    end
    hAnn = annotation('textbox', [0 0.96 1 0.035], ... % 归一化坐标 [x y w h]
        'String','Cluster quality vs K', ...
        'FontWeight','bold','FontSize',12, ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'EdgeColor','none','Interpreter','none');

    exportgraphics(gcf, fullfile(outdir,'ext_quality_vs_K.pdf'), 'ContentType','vector');


    print(fullfile(outdir,'ext_quality_vs_K.png'),'-dpng','-r300');
    print(fullfile(outdir,'ext_quality_vs_K.pdf'),'-dpdf','-r300');
end

%% (c) Cluster temporal prevalence（如有日期）
edges = []; share = [];
if ~isempty(R.timepoint) && ~isempty(R.labels)
    t = R.timepoint(:); L = R.labels(:); K = max(L);
    try
        if ~isdatetime(t), t = datetime(t,'ConvertFrom','datenum'); end
    end
    edges = unique(dateshift(t,'start','month'));
    if numel(edges)<3, edges = unique(dateshift(t,'start','quarter')); end
    if numel(edges)>=3
        share = zeros(numel(edges)-1, K);
        for i=1:numel(edges)-1
            idx = t>=edges(i) & t<edges(i+1); n = sum(idx);
            if n>0, for k=1:K, share(i,k)=sum(L(idx)==k)/n; end, end
        end
        figure('Color','w','Position',[100 100 900 500]);
        area(edges(1:end-1), share); datetick('x','yyyy-mm','keeplimits');
        ylim([0,1]); ylabel('Share'); grid on;
        title('Temporal prevalence of clusters');
        legend(arrayfun(@(k)sprintf('Cluster %d',k),1:K,'UniformOutput',false), ...
            'Location','bestoutside');
        print(fullfile(outdir,'ext_cluster_temporal_prevalence.png'),'-dpng','-r300');
        print(fullfile(outdir,'ext_cluster_temporal_prevalence.pdf'),'-dpdf','-r300');
    end
end

disp('Done. 输出已保存至 figs_ext/ 目录。');

%% ===================== 计算 quality 的子函数 =====================
function [Klist, sil, silU, gap] = calc_quality_byK(S, Win, B)
% 依据 S.byK 循环计算 silhouette / silhouette_unique / gap
assert(isfield(S,'byK') && ~isempty(S.byK), 'S.byK 不存在，无法计算质量曲线');
Klist = sort(arrayfun(@(z) z.K, S.byK));
nK    = numel(Klist);

% 距离矩阵：若 S.D 存在且维度匹配 N，则复用
if isfield(S,'D') && ~isempty(S.D)
    Dfull = S.D;
else
    Rtmp = resolve_result_struct_ext(S);
    fprintf('[Info] 重新计算 DTW 距离（Win=%.2f）…\n', Win);
    Dfull = pairwise_dtw(Rtmp.F_norm, 'Window', Win, 'UseParallel', false);
end
N  = size(Dfull,1);

sil  = nan(nK,1);   % 正常 silhouette
silU = nan(nK,1);   % 去除同省份重复（示例）
gap  = nan(nK,1);   % GAP

for i = 1:nK
    Ki   = Klist(i);
    lab  = S.byK(find([S.byK.K]==Ki,1)).labels(:);
    assert(numel(lab)==N, 'labels 长度与 D 不一致 (K=%d)', Ki);

    sAll = silhouette_precomputed(Dfull, lab);
    sil(i) = mean(sAll, 'omitnan');

    % silhouette_unique：示例按 “省份” 去重 —— 若 meta.province 缺失，则与 sil 相同
    prov = pick(struct_pick(S,'meta'), ["province","prov","Province"]);
    if ~isempty(prov) && numel(prov)==N
        prov = string(prov(:));
        keep = true(N,1);
        for k = 1:Ki
            idx = find(lab==k);
            [~, ia] = unique(prov(idx),'stable');
            keep(idx(setdiff(1:numel(idx), ia))) = false;
        end
        silU(i) = mean(sAll(keep), 'omitnan');
    else
        silU(i) = sil(i);
    end

    % GAP（简易）：对 label 重排 B 次，计算离散度差
    Wk   = mean_intra(Dfull, lab);
    Wkb  = zeros(B,1);
    for b = 1:B
        labB = lab(randperm(N));
        Wkb(b) = mean_intra(Dfull, labB);
    end
    gap(i) = log(Wk) - mean(log(Wkb));
end
end

function W = mean_intra(D, lab)
% 平均簇内距离（upper-tri）
K = max(lab); acc = 0; cnt = 0;
for k = 1:K
    idx = find(lab==k);
    if numel(idx)>=2
        Dk = D(idx,idx); upper = Dk(triu(true(size(Dk)),1));
        acc = acc + sum(upper,'omitnan'); cnt = cnt + numel(upper(~isnan(upper)));
    end
end
W = acc / max(cnt,1);
end

%% ======================== 工具函数（与前版相同） ========================
function R = resolve_result_struct_ext(S)
R = struct;
flds = fieldnames(S);
if numel(flds)==1 && isstruct(S.(flds{1})), S = S.(flds{1}); end
R.expl_var = pick(S,["expl_var","vecr","explained","explained_variance"]);
QQ = struct_pick(S,"cluster_quality"); if isempty(QQ), QQ = struct_pick(S,"quality"); end
R.K_grid            = pick(QQ,["K_grid","Klist","K"]);
R.silhouette        = pick(QQ,["silhouette","silhouette_mean","sil"]);
R.silhouette_unique = pick(QQ,["silhouette_unique","sil_unique"]);
R.gap               = pick(QQ,["gap","gap_stat"]);
M = struct_pick(S,"meta"); if isempty(M), M = struct_pick(S,"metadata"); end
R.timepoint = pick(M,["date","timepoint","t0_date","peak_date"]);
R.labels    = pick(struct_pick(S,"cluster"),["labels","idx"]);
R.F_norm    = pick(S,["F_norm","F","curves","X","Y"]);
R.u         = pick(S,["u","grid","time_grid","U"]);
if size(R.F_norm,1) ~= numel(R.u), R.F_norm = R.F_norm.'; end
end

function val = pick(S, names)
val=[];
if isempty(S), return; end
if istable(S)
    for k=names, k=char(k);
        if any(strcmpi(S.Properties.VariableNames,k)), val=S.(k); return; end
    end
elseif isstruct(S)
    f=fieldnames(S);
    for k=names, k=char(k);
        hit=f(strcmpi(f,k)); if ~isempty(hit), val=S.(hit{1}); return; end
    end
end
end
function S2 = struct_pick(S,k), S2=[]; if isstruct(S)&&isfield(S,k), S2=S.(k); end, end
function ensure_dir(d), if ~exist(d,'dir'), mkdir(d); end, end
function s = silhouette_precomputed(D, labels)
labels = double(labels(:)); N = numel(labels); K = max(labels);
s = nan(N,1);
for i=1:N
    li = labels(i); Di = D(i,:);
    ai = mean(Di(labels==li & (1:N)'~=i),'omitnan');
    bi = inf;
    for k=1:K
        if k==li, continue, end
        bk = mean(Di(labels==k),'omitnan');
        if bk < bi, bi = bk; end
    end
    s(i) = (bi - ai) / max(ai, bi);
end
end
