function create_fig_pairs_panel(scores_root, features_mat, varargin)
% CREATE_FIG_PAIRS_PANEL
% 面板：目标指标 vs 协变量（带分组），含最底部一行组柱状图。
%
% 例：
% create_fig_pairs_panel('scores','features_aug.mat', ...
%   'target','RRI_VA', 'use_rank',true, 'auto_top_k',6, ...
%   'export_pdf',false,'show_fig',true, ...
%   'bar_metric','n','bar_k',3);

% ---------- 参数 ----------
p = inputParser; p.FunctionName = mfilename;
addRequired(p,'scores_root',@(s)ischar(s)||isstring(s));
addRequired(p,'features_mat',@(s)ischar(s)||isstring(s));
addParameter(p,'target','SRS',@(s) any(strcmpi(string(s),{'RRI_VA','MEI_VA','SRS','LPI'})));
addParameter(p,'covars',{},@(c) iscellstr(c) || all(arrayfun(@(z)isstring(z),c))); % 手动指定（白名单中文标签）
addParameter(p,'auto_top_k',6,@(x)isnumeric(x)&&isscalar(x)&&x>=1);               % 自动 Top-K（当 covars 为空时）
addParameter(p,'use_rank',true,@(x)islogical(x)||ismember(x,[0 1]));              % 指标是否先秩化
addParameter(p,'winsor',[0.01 0.99],@(v)isnumeric(v)&&numel(v)==2&&v(1)>=0&&v(2)<=1&&v(1)<v(2));
addParameter(p,'group_mode','region',@(s) any(strcmpi(string(s),{'auto','region','cluster'})));
addParameter(p,'group_k',3,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
addParameter(p,'scatter_alpha',0.65,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<=1);
addParameter(p,'export_pdf',false,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'out_pdf','Fig_Pairs.pdf',@(s)ischar(s)||isstring(s));
addParameter(p,'show_fig',true,@(x)islogical(x)||ismember(x,[0 1]));
% 新增：底部柱状图
addParameter(p,'bar_metric','n',@(s) any(strcmpi(string(s),{'n','median_target','mean_target'})));
addParameter(p,'bar_k',3,@(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>=1)); % []=全部；默认3=最多的三组

addParameter(p,'mark_log',true,@(x)islogical(x)||ismember(x,[0 1]));   % 在散点回归图上标记 log
addParameter(p,'show_shift',true,@(x)islogical(x)||ismember(x,[0 1])); % 平移量是否显示

parse(p,scores_root,features_mat,varargin{:}); A=p.Results;

% ---------- 全局常量/样式 ----------
ABS_TOP3  = true;   % 对角/上三角只取组内 |rho| 最大的前三
MIN_G_N   = 3;      % 计算组内 rho/密度的最小样本阈值
axisCol   = [0.10 0.10 0.10];    % 坐标轴颜色
tickFS    = 10;                  % 刻度字号
textCol   = [55 65 81]/255;      % 文本深灰
COLBASE   = [0 191 166; 58 134 255; 108 76 207]/255; % 对角三条固定基色（青/蓝/紫）

% ---------- 读指标表 ----------
S1 = readtable(fullfile(A.scores_root,'map_ready_SRS.csv'),'VariableNamingRule','preserve');
S2 = readtable(fullfile(A.scores_root,'map_ready_SRS_LPI.csv'),'VariableNamingRule','preserve');
S1.province = string(S1.province); S2.province = string(S2.province);
K = innerjoin(S1(:,{'province','RRI_VA','MEI_VA','SRS'}), S2(:,{'province','LPI'}),'Keys','province');
K = sortrows(K,'province');
if A.use_rank
    for v = ["RRI_VA","MEI_VA","SRS","LPI"], K.(v)=tiedrank(double(K.(v))); end
end
target = upper(string(A.target));

% ---------- 读 features / 分组 ----------
Feat = load_features_aug(A.features_mat);
if ~isstring(Feat.province), Feat.province = string(Feat.province); end
Feat = Feat(:, setdiff(Feat.Properties.VariableNames, {'wave_id'}));
gcol = find_group_col(Feat);
RegMap = unique(Feat(:, ["province", gcol]));
RegMap.Properties.VariableNames{2} = 'group';
RegMap.group = string(RegMap.group);

% 数值列聚合到省级
V = Feat.Properties.VariableNames;
isNum = varfun(@(x) isnumeric(x)||islogical(x), Feat,'OutputFormat','uniform');
isNum(strcmp(V,'province')) = false;
numVarsAll = V(isNum);
G = findgroups(Feat.province);
Agg = table(unique(string(Feat.province)),'VariableNames',{'province'});
for i=1:numel(numVarsAll)
    x = Feat.(numVarsAll{i});
    Agg.(numVarsAll{i}) = splitapply(@(z) median(double(z),'omitnan'), x, G);
end
Agg.province = string(Agg.province);
Agg = sortrows(Agg,'province');

% ---------- 白名单映射 & 协变量选择 ----------
WL  = covar_whitelist_table();
Map = map_covar_to_columns(Agg.Properties.VariableNames, WL); % label/col/is_numeric
numMap = Map([Map.is_numeric]);
labels_all = string({numMap.label}); cols_all = string({numMap.col});

labels = string(A.covars);
if isempty(labels)
    % 优先从 varscreen_<target>.csv 读 Top-K；否则即时计算
    f = fullfile(A.scores_root, sprintf('varscreen_%s.csv', char(target)));
    if isfile(f)
        Tsel = readtable(f,'VariableNamingRule','preserve');
        labels = string(Tsel.label(1:min(A.auto_top_k,height(Tsel))));
    else
        AggW = Agg;
        for k=1:numel(cols_all), v = cols_all(k); AggW.(v) = winsorize_vec(double(AggW.(v)), A.winsor(1), A.winsor(2)); end
        Xtmp = innerjoin(K(:,["province",char(target)]), AggW, 'Keys','province');
        rhoTmp = nan(numel(cols_all),1);
        for k=1:numel(cols_all)
            x = double(Xtmp.(cols_all(k))); y = double(Xtmp.(char(target)));
            m = isfinite(x)&isfinite(y); if sum(m)>=5, rhoTmp(k)=corr(x(m),y(m),'Type','Spearman'); end
        end
        [~,ord]=sort(abs(rhoTmp),'descend'); labels = labels_all(ord(1:min(A.auto_top_k,numel(ord))));
    end
end
colsSel = strings(size(labels));
for i=1:numel(labels)
    idx = find(strcmpi(labels_all, labels(i)),1);
    if isempty(idx), error('协变量“%s”不在白名单映射中。', labels(i)); end
    colsSel(i) = cols_all(idx);
end
for k=1:numel(colsSel)
    v = colsSel(k); Agg.(v) = winsorize_vec(double(Agg.(v)), A.winsor(1), A.winsor(2));
end

% ---------- 合并构造 X（第一列=目标） ----------
T = innerjoin(K, Agg, 'Keys','province');
T = innerjoin(T, RegMap, 'Keys','province'); T.group = string(T.group);
vars = [target; labels(:)];
cols_for_X = [string(target); colsSel(:)];
L = numel(vars);
X = nan(height(T), L);
for j=1:L
    if j==1, X(:,j) = double(T.(char(target))); else, X(:,j) = double(T.(cols_for_X(j))); end
end

% ---------- 分组标签/颜色 ----------
[gLabels, gNames] = build_groups(T, X, A, gcol);
gLabels = gLabels(:);                                      % 统一列向量
Kgrp = numel(gNames);
COLS = palette_colors(Kgrp);                               % 全组配色（用于散点/箱线）
lineCol = COLBASE(min(2,size(COLBASE,1)),:);               % 回归线颜色（基色第2条）

% ---------- 为“每行变量”选择前三强组（优先 |rho|；不足补样本多的组） ----------
top3_by_row = cell(L,1);
grpN_total = accumarray(gLabels,1,[Kgrp,1],@sum,0);
for i=1:L
    y = X(:,1);
    x = X(:,i);
    n_xy = accumarray(gLabels, isfinite(x) & isfinite(y), [Kgrp,1], @sum, 0);
    n_x  = accumarray(gLabels, isfinite(x), [Kgrp,1], @sum, 0);

    rg = nan(Kgrp,1);
    for g=1:Kgrp
        mg = (gLabels==g) & isfinite(x) & isfinite(y);
        if sum(mg) >= MIN_G_N
            rg(g) = corr(x(mg), y(mg), 'Type','Spearman', 'Rows','complete');
        end
    end
    has_rho = isfinite(rg);
    score   = abs(rg); score(~has_rho) = -Inf;
    [~,ord] = sort(score, 'descend');
    pick = ord(score(ord) > -Inf);
    pick = pick(1:min(3,numel(pick)));

    if numel(pick)<3
        [~,ordx] = sort(n_x,'descend');
        fill = setdiff(ordx, pick, 'stable');
        pick = [pick(:).' , fill(1:min(3-numel(pick), numel(fill)))];
    end
    pick = unique(pick,'stable');
    top3_by_row{i} = pick(1:min(3,numel(pick)));
end

% ---------- 整体 rho / p（做 FDR） 与组内 rho ----------
Pair = struct(); P_all=[]; ij_all=[];
for i=1:L
    for j=1:L
        xi=X(:,j); yi=X(:,i); m=isfinite(xi)&isfinite(yi);
        if sum(m)>=5, [R,P]=corr(xi(m),yi(m),'Type','Spearman','Rows','complete'); else, R=NaN; P=NaN; end
        Pair(i,j).rho=R; Pair(i,j).p=P;
        Pair(i,j).rho_g = group_rhos(xi, yi, gLabels, Kgrp, MIN_G_N);
        if i<j && isfinite(P), P_all(end+1,1)=P; ij_all(end+1,:)=[i j]; end %#ok<AGROW>
    end
end
q_all = fdr_bh(P_all);
Qmap = containers.Map("KeyType","char","ValueType","double");
for k=1:numel(P_all), Qmap(sprintf('%d_%d',ij_all(k,1),ij_all(k,2))) = q_all(k); end

% ---------- 出图：多一行给“组柱状图” ----------
nRows = L;  % 底部额外一行
fig = figure('Color','w','Position',[40 40 1280 1020],'Visible','on');
movegui(fig,'center'); plotedit(fig,'on');
tl  = tiledlayout(nRows, L+1, 'Padding','compact', 'TileSpacing','compact');

% 标题（annotation，可拖拽/编辑）
str = sprintf('%s vs covariates (group=%s)', char(target), char(A.group_mode));
annotation(gcf,'textbox',[0.5 0.985 0 0], 'String',str, 'FitBoxToText','on', ...
    'Units','normalized','HorizontalAlignment','center','VerticalAlignment','top', ...
    'LineStyle','none','Interpreter','none','FontWeight','bold',FontSize=15);

% —— 主面板 —— 
for i=1:L
    topG = unique(top3_by_row{i},'stable');
    for j=1:(L+1)
        ax = nexttile(tl, (i-1)*(L+1)+j); hold(ax,'on'); box(ax,'on');
        if j<=L
            if i==j
                % 对角：三条核密度（颜色按名次固定）
                xi = X(:,j); xi = xi(:);
                for k=1:numel(topG)
                    g = topG(k);
                    mg = isfinite(xi) & (gLabels == g);
                    if sum(mg) >= 3
                        [f,xx] = ksdensity(xi(mg));
                        % clr = COLBASE(min(k,size(COLBASE,1)),:);
                        clr = [0.35 0.68 0.87; 0.56 0.81 0.94; 0.00 0.50 0.94];
                        patch(xx, f, clr(k,:), 'FaceAlpha',0.30, ...
                              'EdgeColor',clr(k,:), 'LineWidth',1.0);
                    end
                end
                xl = [min(X(:,j),[],'omitnan'), max(X(:,j),[],'omitnan')];
                if any(~isfinite(xl)), xl=[0 1]; end
                xlim(ax, xl);
                set(ax,'YTick',[], 'XColor',axisCol, 'YColor',axisCol, ...
                       'LineWidth',1.0, 'FontSize',tickFS);
                title(ax, strrep(vars(j),'_','\_'),'Color',textCol,'FontWeight','bold');

            elseif i>j
                % 下三角：散点 + log–log 回归线
                x = X(:,j); y = X(:,i); m = isfinite(x)&isfinite(y);
                [xp, xshift] = to_log_safe(x);  [yp, yshift] = to_log_safe(y);
                mm = m & isfinite(xp) & isfinite(yp);
                for g=1:Kgrp
                    mg = (gLabels==g) & mm;
                    if any(mg)
                        scatter(ax, xp(mg), yp(mg), 16, 'filled', ...
                            'MarkerFaceColor',COLS(g,:), 'MarkerEdgeColor','none', ...
                            'MarkerFaceAlpha',A.scatter_alpha);
                    end
                end
                if sum(mm)>=8
                    try
                        xx = logspace(log10(min(xp(mm))), log10(max(xp(mm))), 100);
                        b  = polyfit(log10(xp(mm)), log10(yp(mm)), 1);
                        yy = 10.^(b(1)*log10(xx) + b(2));
                        plot(ax, xx, yy, '-', 'Color', [0 191 166]/255, 'LineWidth', 2.2);
                    catch, end
                end
                % 仅在最左列/底行标注轴标题，避免拥挤
                if j==1
                    ylabel(ax, sprintf('%s (log_{10})', strrep(vars(i),'_','\_')));
                end
                if i==L && j<=L
                    xlabel(ax, sprintf('%s (log_{10})', strrep(vars(j),'_','\_')));
                end

                set(ax,'XScale','log','YScale','log', ...
                    'XColor',axisCol,'YColor',axisCol, ...
                    'LineWidth',1.0,'FontSize',tickFS);
                % —— 打上 log 变换标记（可选显示平移量） ——
                if A.mark_log
                    tag = 'log_{10}–log_{10}';
                    if A.show_shift && (xshift>0 || yshift>0)
                        % 显示为：log10(x + a), log10(y + b)
                        tag = sprintf('%s (x+%g, y+%g)', tag, xshift, yshift);
                    end
                    text(ax, 0.03, 0.08, tag, 'Units','normalized', ...
                        'Color',[0.30 0.30 0.30], 'FontSize',8, 'Interpreter','tex');
                end


            else
                % 上三角：整体ρ + FDR(q) + 前三组ρ
                % —— 上三角：整体ρ + FDR(q) + 前三组ρ（固定三位小数并带 n） ——
                R  = Pair(i,j).rho; P = Pair(i,j).p; rg = Pair(i,j).rho_g;

                % 计算这对在各组的“共同有效样本量”（x、y 同时非缺失）
                xi = X(:,j); yi = X(:,i);
                n_xy = accumarray(gLabels, isfinite(xi) & isfinite(yi), [Kgrp,1], @sum, 0);

                key = sprintf('%d_%d', i, j);
                qv = NaN; if isKey(Qmap,key), qv = Qmap(key); end

                % 固定三位小数的格式
                fmtR = '%.3f';  fmtQ = '%.3g';  % q 用工程记数更友好

                txt = sprintf(['Corr: ' fmtR '%s\nFDR q: ' fmtQ], R, stars(P), qv);

                % 只列出“该行挑中的前三组”
                topG = top3_by_row{i};
                for k = 1:numel(topG)
                    g = topG(k);
                    if isfinite(rg(g))
                        txt = sprintf('%s\n%s: %.1f', txt, gNames{g}, rg(g));
                    end
                end

                text(ax, 0.02,0.95, txt, 'Units','normalized','VerticalAlignment','top', ...
                    'Color', textCol,'FontWeight','bold','Interpreter','none');
                axis(ax,'off');

            end
        else
            % 最右列：该行变量的分组箱线（全部组）
            xi = X(:,i); mgAll = isfinite(xi) & isfinite(gLabels);
            if any(mgAll)
                for gi=1:Kgrp
                    mg = (gLabels==gi) & mgAll;
                    if any(mg)
                        bc = boxchart(ax, gi*ones(sum(mg),1), xi(mg));
                        bc.BoxFaceColor = COLS(gi,:); bc.BoxFaceAlpha = 1.0;
                        bc.BoxEdgeColor = axisCol;     bc.MarkerStyle = 'none';
                    end
                end
                gShort = cellfun(@(s) mk_abbrev(s), gNames, 'UniformOutput', false);
                xlim(ax,[0.5, Kgrp+0.5]);
                set(ax,'XTick',1:Kgrp,'XTickLabel',gShort,'XTickLabelRotation',45, ...
                       'XColor',axisCol,'YColor',axisCol,'LineWidth',1.0,'FontSize',tickFS);
            else
                axis(ax,'off');
            end
        end
        if j>1 && j<=L, set(ax,'YTickLabel',[]); end
        if i<L  && j<=L, set(ax,'XTickLabel',[]); end
    end
end

% % —— 最底部一行：组柱状图（跨整行一个格） —— 
% axB = nexttile(tl, (nRows-1)*(L+1)+1, [1, L+1]); hold(axB,'on'); box(axB,'on');
% switch lower(char(A.bar_metric))
%     case 'n'
%         ybar = accumarray(gLabels,1,[Kgrp,1],@sum,0);
%         ylab = 'Number of provinces';
%     case 'median_target'
%         ybar = splitapply(@(z) median(z,'omitnan'), X(:,1), gLabels);
%         ylab = sprintf('Median of %s', char(target));
%     case 'mean_target'
%         ybar = splitapply(@(z) mean(z,'omitnan'), X(:,1), gLabels);
%         ylab = sprintf('Mean of %s', char(target));
%     otherwise
%         ybar = accumarray(gLabels,1,[Kgrp,1],@sum,0);
%         ylab = 'Number of provinces';
% end
% % 若 bar_k 指定仅画前 K（按 ybar 降序）
% ordBar = (1:Kgrp)';
% if ~isempty(A.bar_k)
%     [~,ord] = sort(ybar,'descend');
%     ordBar = ord(1:min(A.bar_k, numel(ord)));
% end
% bcol = COLS(ordBar,:); 
% bar(axB, 1:numel(ordBar), ybar(ordBar), 0.7, 'FaceColor','flat','CData',bcol, 'EdgeColor',axisCol);
% gShort = cellfun(@(s) mk_abbrev(s), gNames, 'UniformOutput', false);
% set(axB,'XTick',1:numel(ordBar),'XTickLabel',gShort(ordBar), ...
%          'XColor',axisCol,'YColor',axisCol,'LineWidth',1.0,'FontSize',tickFS);
% ylabel(axB, ylab, 'Color',textCol);
% title(axB, 'Group summary','Color',textCol,'FontWeight','bold');

% 导出
if A.export_pdf
    set(fig,'Renderer','painters');
    exportgraphics(fig, A.out_pdf, 'ContentType','vector');
    fprintf('[export] saved: %s\n', A.out_pdf);
end
if A.show_fig, figure(fig); drawnow; end
end

% ================= helpers =================

function [xp, shift] = to_log_safe(x)
x = double(x); m = isfinite(x);
if ~any(m), xp=x; shift=0; return; end
xmin = min(x(m));
if xmin > 0, xp = x; shift = 0;
else, shift = max(1e-6, 1 - xmin); xp = x + shift; xp(~isfinite(x)) = NaN;
end
end

function rg = group_rhos(x, y, gLabels, K, MIN_G_N)
if nargin<5, MIN_G_N=3; end
rg = nan(K,1); m0 = isfinite(x) & isfinite(y);
for g=1:K
    mg = (gLabels==g) & m0;
    if sum(mg) >= MIN_G_N
        rg(g) = corr(x(mg), y(mg), 'Type','Spearman', 'Rows','complete');
    end
end
end

function q = fdr_bh(p)
q = nan(size(p)); ix=find(isfinite(p)); if isempty(ix), return; end
[sp,ord]=sort(p(ix)); m=numel(sp);
qv = sp .* m ./ (1:m)'; qv = flipud(cummin(flipud(qv)));
q(ix(ord)) = qv;
end

function s = stars(p)
if ~isfinite(p), s=''; return; end
if p<1e-3, s='***'; elseif p<1e-2, s='**'; elseif p<0.05, s='*'; else, s=''; end
end

function xw = winsorize_vec(x, plo, phi)
x = double(x); qlo = quantile(x,plo,'all'); qhi = quantile(x,phi,'all');
xw = x; xw(x<qlo)=qlo; xw(x>qhi)=qhi;
end

function gcol = find_group_col(Feat)
V = string(Feat.Properties.VariableNames); L = lower(V);
if any(contains(L,["regioncategory","region_cat","macroregion","regionclass"]))
    gcol = V(find(contains(L,["regioncategory","region_cat","macroregion","regionclass"]),1,'first')); return;
end
if any(contains(L,["climatezone","climate_zone","koppen","climate_cat"]))
    gcol = V(find(contains(L,["climatezone","climate_zone","koppen","climate_cat"]),1,'first')); return;
end
if any(contains(L,"region"))
    gcol = V(find(contains(L,"region"),1,'first')); return;
end
error('features_aug 中未找到 region/climate/region category 任一分组列');
end

function [gLabels, gNames] = build_groups(T, X, A, gcol)
switch lower(char(A.group_mode))
    case {'auto','region'}
        g = string(T.group); u = unique(g,'stable'); u(u=="") = [];
        if isempty(u), gLabels = ones(height(T),1); gNames={'All'};
        else, gNames = cellstr(u); gLabels = zeros(height(T),1); for i=1:numel(u), gLabels(g==u(i))=i; end
        end
    case 'cluster'
        Z=X; for j=1:size(Z,2), Z(:,j)=(Z(:,j)-mean(Z(:,j),'omitnan'))/std(Z(:,j),'omitnan'); end
        m=all(isfinite(Z),2); Zc=Z(m,:); [idx,~]=kmeans(Zc, A.group_k,'Replicates',10,'Distance','sqeuclidean');
        gLabels=zeros(height(Z),1); gLabels(m)=idx; gLabels(~m)=mode(gLabels(gLabels>0));
        gNames=arrayfun(@(k)sprintf('Cluster%d',k),1:A.group_k,'UniformOutput',false);
    otherwise, error('group_mode 仅支持 auto/region/cluster');
end
end

function COL = palette_colors(K)
% 全组配色：青绿→蓝→紫插值
base = [0 191 166; 58 134 255; 108 76 207]/255;
if K<=size(base,1), COL = base(1:K,:);
else, COL = interp1(1:3, base, linspace(1,3,K),'linear');
end
end

function WL = covar_whitelist_table()
WL = struct('label',{},'aliases',{},'is_numeric',{});
add("age65 share",                ["age65share","age65","elderlyshare","share65","pop65","age_65"], true);
add("education share",            ["educationshare","edu_share","tertiaryedu","college_share","higheredu"], true);
add("gdp_{pc}",                   ["gdppc","gdppercapita","gdp_pc","log_gdp_pc","gdp_percap"], true);
add("hospital beds per 1k",       ["hospitalbeds","bedsper1000","beds_per1k","beds_pc","bed_density"], true);
add("proportion of students",     ["proportion of students","studentshare","proportionstudents","enrollment_share"], true);
add("inbound overnight visitors", ["overnightvisitors","inboundvisitors","tourist_night","visitor_night"], true);
add("disposable income_{pc}",     ["dispincpc","incomepercapita","disposable_income_pc","income_pc","log_income_pc"], true);
add("passenger volume",           ["passengervolume","passenger_flow","air_passenger","rail_passenger"], true);
add("public transit volume",      ["transitvolume","publictransit","log_transit_vol","metro_bus_rail"], true);
add("region category",            ["regioncategory","region_cat","macroregion","regionclass"], false);
add("tertiary-sector share",      ["tertiaryshare","serviceshare","service_share","tertiary"], true);
add("urban population density",   ["urbanpopulationdensity","urban_density","pop_density_urban"], true);
add("urban rate",                 ["urbanrate","urbanization","urban_ratio","urban_share"], true);
add("climate zone",               ["climatezone","climate_zone","koppen","climate_cat"], false);
    function add(lbl, als, isn), WL(end+1)=struct('label',lbl,'aliases',als,'is_numeric',logical(isn)); end
end

function Map = map_covar_to_columns(colNames, WL)
% 把白名单标签映射到 features 列名（最相近匹配，容错）
C = string(colNames); Ccanon = lower(regexprep(C,'[^a-z0-9]',''));
Map = repmat(struct('label',"",'col',"",'is_numeric',true),0,1);
for i=1:numel(WL)
    lab=string(WL(i).label);
    if ~WL(i).is_numeric
        Map(end+1)=struct('label',lab,'col',lab,'is_numeric',false); %#ok<AGROW>
        continue
    end
    aliases=string(WL(i).aliases); best=""; scoreBest=-inf;
    for j=1:numel(C)
        sc=-inf;
        for a=aliases
            tok=string(regexp(lower(regexprep(a,'[^a-z0-9]','')),'[a-z0-9]+','match'));
            if all(contains(Ccanon(j), tok)), sc=max(sc, numel(tok)+0.001*strlength(a)); end
        end
        if sc>scoreBest, scoreBest=sc; best=C(j); end
    end
    if scoreBest>0
        Map(end+1)=struct('label',lab,'col',best,'is_numeric',true); %#ok<AGROW>
    else
        warning('白名单变量[%s]未在 features 中匹配到列名。', lab);
    end
end
% 去重（用字符串键避免维度不兼容）
cols = string({Map.col})';
isnum = string(double([Map.is_numeric]'));
keys = cols + "__" + isnum;
[~, ia] = unique(keys, 'stable'); Map = Map(ia);
end

function ab = mk_abbrev(name)
% "East China" -> "E.C."；"Southwest Region" -> "S.R."
parts = regexp(name, '[A-Za-z]+', 'match');
if isempty(parts), ab = name; return; end
ini = upper(cellfun(@(p)p(1), parts, 'UniformOutput', false));
ab = strjoin( cellfun(@(c)[c '.'], ini, 'UniformOutput', false), '' );
end
