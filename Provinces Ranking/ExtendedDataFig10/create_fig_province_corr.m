function create_fig_province_corr(scores_root, varargin)
% 省×省相关性热图（SRS/RRI/MEI/LPI）
% - 三角半幅位置由 corner 控制；默认 'SW' = 左下角（西南角）
% - 色标在底部；按层次聚类画三角黑框
% - 可选在图上写每簇指标：r_within、r_between、silhouette
%
% 例：
% create_fig_province_corr('scores','metric','RRI','method','spearman', ...
%     'corner','SW','k',[], 'show_metrics',true, 'export_pdf',true);

p = inputParser; p.FunctionName = mfilename;
addRequired(p,'scores_root',@(s)ischar(s)||isstring(s));
addParameter(p,'metric','SRS',@(s) any(strcmpi(string(s),{'SRS','RRI','MEI','LPI'})));
addParameter(p,'method','spearman',@(s) any(strcmpi(string(s),{'spearman','pearson'})));
addParameter(p,'corner','NW',@(s) any(strcmpi(string(s),{'SW','SE','NW','NE'}))); % 三角放置角
addParameter(p,'lpi_boot_path','scores',@(s)ischar(s)||isstring(s));
addParameter(p,'min_overlap',2,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
addParameter(p,'k',[],@(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>=2)); % 簇数；空=自动阈值
addParameter(p,'show_metrics',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'export_pdf',false,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'out_pdf','Fig_ProvCorr.pdf',@(s)ischar(s)||isstring(s));
parse(p,scores_root,varargin{:}); A=p.Results;

% ---------- 读入并构造省×样本矩阵 ----------
m = upper(string(A.metric));
switch m
    case "SRS"
        T = readtable(fullfile(A.scores_root,'province_trajectories.csv'), 'VariableNamingRule','preserve');
        needCols = {'province','order','SRS_wave'};
        miss = setdiff(needCols, T.Properties.VariableNames);
        if ~isempty(miss), error('province_trajectories.csv 缺少列：%s', strjoin(miss, ', ')); end
        T = T(:, needCols); T.province = string(T.province);
        T.order = normalize_order(T.province, T.order);
        T = sortrows(T,{'province','order'});
        wide = unstack(T, 'SRS_wave', 'order', 'AggregationFunction',@mean, ...
                       'VariableNamingRule','preserve');
        prov = string(wide.province); X = wide{:,2:end};

    case {"RRI","MEI"}
        W0 = readtable(fullfile(A.scores_root,'wave_scores.csv'), 'VariableNamingRule','preserve');
        W0.province = string(W0.province); W0.wave_id = string(W0.wave_id);
        W0.order = normalize_order(W0.province, W0.wave_id);
        ycol = 'RRI_wave'; if strcmpi(m,'MEI'), ycol='MEI_wave'; end
        Wmini = W0(:,{'province','order'}); Wmini.(ycol) = W0.(ycol);
        Wmini = sortrows(Wmini,{'province','order'});
        WIDE  = unstack(Wmini, ycol, 'order', 'AggregationFunction',@mean, 'VariableNamingRule','preserve');
        [prov, X] = wide_to_matrix(WIDE);

    case "LPI"
        path = char(A.lpi_boot_path);
        if isempty(path)
            for c = ["lpi_boot.mat","lpi_boot.csv","LPI_boot.mat","LPI_boot.csv"]
                f = fullfile(A.scores_root, c); if isfile(f), path = f; break; end
            end
        end
        if isempty(path), error('LPI 需要 P×B 样本，请提供 lpi_boot.(mat|csv)。'); end
        [prov, X] = load_boot_matrix(path);

    otherwise
        error('Unknown metric: %s', m);
end

% ---------- 相关矩阵 ----------
mu = nanmean(X,2); sd = nanstd(X,0,2); sd(~isfinite(sd)|sd==0) = 1;
Xz = (X - mu) ./ sd;
P = numel(prov); R = nan(P);
for i=1:P
    for j=i:P
        msk = isfinite(Xz(i,:)) & isfinite(Xz(j,:));
        if sum(msk) >= A.min_overlap
            R(i,j) = corr(Xz(i,msk)', Xz(j,msk)', 'Type',lower(A.method));
            R(j,i) = R(i,j);
        end
    end
end

% ---------- 层次聚类 & 排序 ----------
Rfill = R; Rfill(isnan(Rfill)) = 0;
D = 1 - Rfill; D(eye(P)==1) = 0;
Z = linkage(squareform(max(D,0)), 'average');      % average linkage
ord = optimalleaforder(Z, D);
Rsorted = R(ord,ord); prov_sorted = prov(ord);

% ---------- 按 corner 取半幅 & 设置轴向 ----------
[maskTri, ydir] = triangle_mask_and_dir(P, upper(string(A.corner)));
Rtri = Rsorted; Rtri(~maskTri) = NaN;

% ---------- 分簇（用于黑框 + 指标） ----------
if isempty(A.k)
    cutoff = 0.7 * max(Z(:,3));
    grp = cluster(Z, 'cutoff', cutoff, 'criterion','distance');
else
    grp = cluster(Z, 'maxclust', A.k);
end
grp = grp(ord);

% ---------- 画图（色标在底部） ----------
fig = figure('Color','w','Position',[80 60 1000 900]); movegui(fig,'center'); plotedit(fig,'on');
ax  = axes('Parent',fig,'Position',[0.12 0.12 0.76 0.78]); hold(ax,'on');

hImg = imagesc(ax, Rtri, [-1 1]);                     % 如需更接近参考图，可改 [-0.2 1]
set(hImg,'AlphaData', double(~isnan(Rtri)));
set(ax,'YDir', ydir);                                 % 'normal' 或 'reverse'
colormap(ax, diverging_bluered_local(255));

% 对角线
plot(ax, 0.5+(1:P), 0.5+(1:P), '-', 'Color',[0 0 0], 'LineWidth',1.2);

% 轴与标签
set(ax,'XLim',[0.5 P+0.5], 'YLim',[0.5 P+0.5], 'Box','off','Layer','top','TickLength',[0 0]);
set(ax,'XTick',1:P,'XTickLabel',prov_sorted,'XTickLabelRotation',90, ...
        'YTick',1:P,'YTickLabel',prov_sorted);
ax.LooseInset = max(ax.TightInset, 0.02);

% 底部色标
cb = colorbar(ax,'Location','southoutside'); cb.Label.String = 'Correlation \rho';

% 聚类三角黑框（按 corner 画三条边）
draw_cluster_triangles_corner(ax, grp, upper(string(A.corner)));

% 指标写在图上（可关）
if A.show_metrics
    S = cluster_stats(Rsorted, grp);   % cluster, n, r_within, r_between, sil_mean
    annotate_cluster_stats(ax, S);
end

% 导出
if A.export_pdf
    set(fig,'Renderer','painters');
    exportgraphics(fig, A.out_pdf, 'ContentType','vector');
    fprintf('[export] saved: %s\n', A.out_pdf);
end
end

% ===== helpers =====
function ord = normalize_order(prov, orderCol)
if isnumeric(orderCol), ord = orderCol; return; end
ordNum = str2double(string(orderCol));
if any(~isnan(ordNum)), ord = ordNum;
else, p = string(prov); ord = zeros(size(p)); u = unique(p,'stable');
     for i=1:numel(u), idx=find(p==u(i)); ord(idx)=1:numel(idx); end
end
end

function [prov, X] = wide_to_matrix(W)
prov = string(W.province);
isNum = varfun(@(v) isnumeric(v)||islogical(v), W, 'OutputFormat','uniform'); isNum(1)=false;
if ~any(isNum), error('pivot 结果未产生任何数值列，请检查 order 是否有效。'); end
X = W{:, isNum};
end

function [prov, X] = load_boot_matrix(f)
[~,~,ext]=fileparts(f);
switch lower(ext)
    case '.mat'
        S = load(f); prov=[]; X=[];
        if isfield(S,'lpi_boot') && isnumeric(S.lpi_boot), X=S.lpi_boot; prov="P"+string(1:size(X,1)); end
        if isempty(X) && isfield(S,'Boot') && isfield(S.Boot,'lpi_boot')
            X = S.Boot.lpi_boot; if isfield(S.Boot,'province_order'), prov=string(S.Boot.province_order); else, prov="P"+string(1:size(X,1)); end
        end
        if isempty(X), error('MAT 文件中未找到 P×B 的 LPI 样本。'); end
    otherwise
        T = readtable(f,'VariableNamingRule','preserve');
        assert(ismember('province', T.Properties.VariableNames),'LPI CSV 需含列 province');
        prov = string(T.province);
        X = T{:, setdiff(T.Properties.VariableNames, {'province'})};
end
end

function C = diverging_bluered_local(N)
if nargin<1, N=255; end
t = linspace(0,1,N)'; blue=[0.16 0.36 0.86]; white=[1 1 1]; red=[0.82 0.17 0.17];
C = zeros(N,3);
for i=1:N
    if t(i)<=0.5, C(i,:)=blue+(white-blue)*(t(i)/0.5);
    else,         C(i,:)=white+(red-white)*((t(i)-0.5)/0.5);
    end
end
end

function [maskTri, ydir] = triangle_mask_and_dir(P, corner)
% 返回半幅遮罩与 YDir，确保三角位于指定角
switch upper(string(corner))
    case "SW"   % 左下角：默认
        maskTri = tril(true(P),0); ydir='normal';
    case "NE"   % 右上角
        maskTri = triu(true(P),0); ydir='normal';
    case "NW"   % 左上角
        maskTri = tril(true(P),0); ydir='reverse';
    case "SE"   % 右下角
        maskTri = triu(true(P),0); ydir='reverse';
    otherwise
        error('corner 取值应为 SW/SE/NW/NE');
end
end

function draw_cluster_triangles_corner(ax, grp, corner)
% 按 corner 在三角热图边勾勒聚类边框（三条边）——修正 YDir 影响
[~,~,gidx] = unique(grp,'stable');
P = numel(gidx); hold(ax,'on');
edges = [1; find(diff(gidx)~=0)+1; P+1];

% 关键：根据轴的 YDir 计算“视觉上的”top/bottom 所对应的 y 值
isRev = strcmpi(get(ax,'YDir'),'reverse');

for k = 1:numel(edges)-1
    a = edges(k); b = edges(k+1)-1;    % 簇覆盖 [a..b]
    x0 = a - 0.5; x1 = b + 0.5;
    y0 = a - 0.5; y1 = b + 0.5;

    % 在视觉坐标中，顶部/底部应是哪一个 y？
    yTop    = tern(isRev, y1, y0);     % 视觉的“上边”
    yBottom = tern(isRev, y0, y1);     % 视觉的“下边”

    switch upper(string(corner))
        case "SW"   % 左下角：左边 + 底边 + 对角
            plot(ax,[x0 x0],[y0 y1],'k-','LineWidth',1.2);         % 左边（不受 YDir 影响）
            plot(ax,[x0 x1],[yBottom yBottom],'k-','LineWidth',1.2);% 底边（受 YDir）
            plot(ax,[x0 x1],[y0 y1],'k-','LineWidth',1.2);          % 对角
        case "SE"   % 右下角：右边 + 底边 + 对角
            plot(ax,[x1 x1],[y0 y1],'k-','LineWidth',1.2);         % 右边
            plot(ax,[x0 x1],[yBottom yBottom],'k-','LineWidth',1.2);% 底边（受 YDir）
            plot(ax,[x0 x1],[y0 y1],'k-','LineWidth',1.2);          % 对角
        case "NW"   % 左上角：左边 + 顶边 + 对角
            plot(ax,[x0 x0],[y0 y1],'k-','LineWidth',1.2);         % 左边
            plot(ax,[x0 x1],[yTop yTop],'k-','LineWidth',1.2);      % 顶边（受 YDir）
            plot(ax,[x0 x1],[y0 y1],'k-','LineWidth',1.2);          % 对角
        case "NE"   % 右上角：右边 + 顶边 + 对角
            plot(ax,[x1 x1],[y0 y1],'k-','LineWidth',1.2);         % 右边
            plot(ax,[x0 x1],[yTop yTop],'k-','LineWidth',1.2);      % 顶边（受 YDir）
            plot(ax,[x0 x1],[y0 y1],'k-','LineWidth',1.2);          % 对角
        otherwise
            error('corner 取值应为 SW/SE/NW/NE');
    end
end
end

function y = tern(cond,a,b), if cond, y=a; else, y=b; end, end


function S = cluster_stats(Rsorted, grp)
% 每簇：数量 n、组内平均相关 r_within、组间平均相关 r_between、平均轮廓系数 sil_mean
[~,~,gidx]=unique(grp,'stable'); P=numel(gidx); K=max(gidx);
S = table((1:K)', zeros(K,1), zeros(K,1), zeros(K,1), zeros(K,1), ...
          'VariableNames', {'cluster','n','r_within','r_between','sil_mean'});
D = 1 - Rsorted; D(eye(P)==1) = 0;

% silhouette（优先用预计算距离）
sil = NaN(P,1);
try
    sil = silhouette(D, gidx, 'Distance','precomputed');
catch
    for i=1:P
        ci = gidx(i); in=(gidx==ci); out=~in;
        ai = mean(D(i, in & (1:P)'~=i), 'omitnan');
        bj = inf; for c=1:K, if c==ci, continue; end, bj=min(bj, mean(D(i,gidx==c),'omitnan')); end
        sil(i) = (bj - ai) / max(ai, bj);
    end
end

for c=1:K
    m = (gidx==c); S.n(c)=sum(m);
    Rc = Rsorted(m,m); Rc(eye(sum(m))==1) = NaN;
    S.r_within(c)  = mean(Rc,'all','omitnan');
    Ro = Rsorted(m,~m);
    S.r_between(c) = mean(Ro,'all','omitnan');
    S.sil_mean(c)  = mean(sil(m),'omitnan');
end
end

function annotate_cluster_stats(ax, S)
% 在图右下角写簇指标
lines = arrayfun(@(i) sprintf('C%d: n=%d, r_w=%.2f, r_b=%.2f, sil=%.2f', ...
                S.cluster(i), S.n(i), S.r_within(i), S.r_between(i), S.sil_mean(i)), ...
                (1:height(S))', 'UniformOutput', false);
txt = strjoin(lines, '\n');
xl = xlim(ax); yl = ylim(ax);
xpos = xl(1) + 0.62*(xl(2)-xl(1));
ypos = yl(1) + 0.10*(yl(2)-yl(1));
t = text(ax, xpos, ypos, txt, 'FontSize',9, 'VerticalAlignment','bottom', ...
         'BackgroundColor',[1 1 1 0.9], 'Margin',6, 'EdgeColor',[0.2 0.2 0.2]);
uistack(t,'top');
end
