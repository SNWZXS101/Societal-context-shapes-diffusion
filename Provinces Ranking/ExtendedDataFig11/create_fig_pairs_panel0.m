function create_fig_pairs_panel0(scores_root, features_mat, varargin)
% CREATE_FIG_PAIRS_PANEL  一张“相关性面板”：【一个指标 + 白名单协变量 Top-K】
% - 若未指定 covars，则自动读取 scores/varscreen_<TARGET>.csv 的 Top-K；
%   若文件不存在，则即时计算 |rho| 选 Top-K（仅在白名单数值项内）。
% - 对角：核密度（分组叠加）；下三角：散点+线性趋势；上三角：相关文本（整体 + 分组ρ）；
% - 右侧一列：按组箱线；分组优先 'region category'→'climate zone'→features里的'region'；
% - 默认 rank 版（Spearman 等价）；配色=日照金山后半段（青绿/蓝/紫）；作图后强制 FIG 弹窗。

% ========= 参数 =========
p = inputParser; p.FunctionName = mfilename;
addRequired(p,'scores_root',@(s)ischar(s)||isstring(s));
addRequired(p,'features_mat',@(s)ischar(s)||isstring(s));
addParameter(p,'target','SRS',@(s) any(strcmpi(string(s),{'RRI_VA','MEI_VA','SRS','LPI'})));
addParameter(p,'covars',{},@(c) iscellstr(c) || all(arrayfun(@(z)isstring(z),c))); % 白名单里的 label
addParameter(p,'auto_top_k',6,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'use_rank',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'winsor',[0.01 0.99],@(v)isnumeric(v)&&numel(v)==2&&v(1)>=0&&v(2)<=1&&v(1)<v(2));
addParameter(p,'group_mode','auto',@(s) any(strcmpi(string(s),{'auto','region','cluster'})));
addParameter(p,'group_k',3,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
addParameter(p,'scatter_alpha',0.68,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<=1);
addParameter(p,'annot_sig','stars',@(s) any(strcmpi(string(s),{'none','stars'})));
addParameter(p,'palette','alpenglow_cool',@(s)ischar(s)||isstring(s));
addParameter(p,'export_pdf',false,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'out_pdf','Fig_Pairs.pdf',@(s)ischar(s)||isstring(s));
addParameter(p,'show_fig',true,@(x)islogical(x)||ismember(x,[0 1]));
parse(p,scores_root,features_mat,varargin{:}); A=p.Results;

% ========= 指标 =========
S1 = readtable(fullfile(A.scores_root,'map_ready_SRS.csv'), 'VariableNamingRule','preserve');
S2 = readtable(fullfile(A.scores_root,'map_ready_SRS_LPI.csv'), 'VariableNamingRule','preserve');
S1.province = string(S1.province); S2.province = string(S2.province);
K = innerjoin(S1(:,{'province','RRI_VA','MEI_VA','SRS'}), S2(:,{'province','LPI'}), 'Keys','province');
K = sortrows(K,'province'); if A.use_rank, for v=["RRI_VA","MEI_VA","SRS","LPI"], K.(v)=tiedrank(double(K.(v))); end, end
target = upper(string(A.target));

axisCol = [0.30 0.30 0.30];  % 深灰轴色（建议放在文件顶部全局复用）
tickFS  = 9;                 % 刻度字号（全局）


% ========= features（聚合） + 分组列 =========
Feat = load_features_aug(A.features_mat);
if ~isstring(Feat.province), Feat.province = string(Feat.province); end
Feat = Feat(:, setdiff(Feat.Properties.VariableNames, {'wave_id'}));
gcol = find_group_col(Feat);
RegMap = unique(Feat(:, ["province", gcol]));   % 统一用 string 索引
RegMap.Properties.VariableNames{2}='group'; RegMap.group=string(RegMap.group);

V = Feat.Properties.VariableNames;
isNum = varfun(@(x) isnumeric(x)||islogical(x), Feat,'OutputFormat','uniform'); isNum(strcmp(V,'province'))=false;
numVarsAll = V(isNum);
G = findgroups(Feat.province);
Agg = table(unique(string(Feat.province)),'VariableNames',{'province'});
for i=1:numel(numVarsAll)
    x = Feat.(numVarsAll{i});
    Agg.(numVarsAll{i}) = splitapply(@(z) median(double(z),'omitnan'), x, G);
end
Agg.province = string(Agg.province);
Agg = sortrows(Agg,'province');

% ========= 白名单映射（仅数值项） =========
WL = covar_whitelist_table();
Map = map_covar_to_columns(Agg.Properties.VariableNames, WL);
numMap = Map([Map.is_numeric]);
labels_all = string({numMap.label}); cols_all = string({numMap.col});

% ========= 协变量选择：若未指定，则从 varscreen/<target> 或即时计算中挑 Top-K =========
labels = string(A.covars);
if isempty(labels)
    f = fullfile(A.scores_root, sprintf('varscreen_%s.csv', char(target)));
    if isfile(f)
        T = readtable(f,'VariableNamingRule','preserve');
        labels = string(T.label(1:min(A.auto_top_k, height(T))));
    else
        % 即时计算 |rho| 选 Top-K（仅白名单）
        AggW = Agg;
        for k=1:numel(cols_all), v = cols_all(k); AggW.(v) = winsorize_vec(double(AggW.(v)), A.winsor(1), A.winsor(2)); end
        X = innerjoin(K(:,["province",char(target)]), AggW, 'Keys','province');
        X = innerjoin(X, RegMap, 'Keys','province'); %#ok<NASGU>
        RHO = nan(numel(cols_all),1);
        for k=1:numel(cols_all)
            x = double(X.(cols_all(k))); y = double(X.(char(target)));
            m = isfinite(x)&isfinite(y); if sum(m)>=5, RHO(k)=corr(x(m),y(m),'Type','Spearman'); end
        end
        [~,ord]=sort(abs(RHO),'descend'); labels = labels_all(ord(1:min(A.auto_top_k, numel(ord))));
    end
end

% 将 labels → 实际列名
cols = strings(size(labels));
for i=1:numel(labels)
    idx = find(strcmpi(labels_all, labels(i)),1);
    if isempty(idx), error('协变量“%s”不在白名单映射中。', labels(i)); end
    cols(i) = cols_all(idx);
end

% ========= Winsorize + 对齐 =========
for k=1:numel(cols)
    v = cols(k); Agg.(v) = winsorize_vec(double(Agg.(v)), A.winsor(1), A.winsor(2));
end
T = innerjoin(K, Agg, 'Keys','province');
T = innerjoin(T, RegMap, 'Keys','province');
T.group = string(T.group);

% ========= 构造绘图矩阵 =========
vars = [target; labels(:)];
cols_for_X = [string(target); cols(:)];
L = numel(vars);
X = nan(height(T), L);
for j=1:L
    if j==1, X(:,j) = double(T.(char(target)));
    else,     X(:,j) = double(T.(cols_for_X(j)));
    end
end

% ========= 分组 =========
[gLabels, gNames] = build_groups(T, X, A, gcol);
col = palette_colors(A.palette, numel(gNames));
fg = [55 65 81]/255;

% ========= 作图（保证 FIG 弹窗） =========
fig = figure('Color','w','Position',[50 50 1200 980],'Visible','on'); movegui(fig,'center'); plotedit(fig,'on');
tl  = tiledlayout(L, L+1, 'Padding','compact', 'TileSpacing','compact');

% 预计算 Spearman（整体 & 分组）
Pair = struct();
for i=1:L
    for j=1:L
        xi=X(:,j); yi=X(:,i); m=isfinite(xi)&isfinite(yi);
        if sum(m)>=5, [R,P]=corr(xi(m),yi(m),'Type','Spearman','Rows','complete'); else, R=NaN; P=NaN; end
        rg = nan(numel(gNames),1);
        for g=1:numel(gNames)
            mg=m & gLabels==g; if sum(mg)>=5, rg(g)=corr(xi(mg),yi(mg),'Type','Spearman'); end
        end
        Pair(i,j).rho=R; Pair(i,j).p=P; Pair(i,j).rho_g=rg;
    end
end

for i=1:L
    for j=1:(L+1)
        nexttile(tl, (i-1)*(L+1)+j); hold on; axis tight; box off;
        if j<=L
            if i==j
                xi=X(:,j);
                for g=1:numel(gNames)
                    mg=isfinite(xi)&gLabels==g;
                    if sum(mg)>=5
                        [f,xx]=ksdensity(xi(mg));
                        patch(xx,f,col(g,:), 'FaceAlpha',0.52,'EdgeColor','none');
                    end
                end
                xl=[min(X(:,j),[],'omitnan'), max(X(:,j),[],'omitnan')]; if any(~isfinite(xl)), xl=[0 1]; end
                xlim(xl); 
                set(gca,'XColor',axisCol,'YColor',axisCol);
                title(strrep(vars(j),'_','\_'),'Color',fg,'FontWeight','bold');
            elseif i>j
                x=X(:,j); yv=X(:,i); m=isfinite(x)&isfinite(yv);
                for g=1:numel(gNames)
                    mg=m & gLabels==g;
                    if any(mg)
                        scatter(x(mg), yv(mg), 12,'filled','MarkerFaceColor',col(g,:), ...
                            'MarkerEdgeColor','none','MarkerFaceAlpha',A.scatter_alpha);
                    end
                end
                if sum(m)>=8
                    try
                        xx=linspace(min(x(m)),max(x(m)),100);
                        pf=polyfit(x(m), yv(m), 1); yline=polyval(pf, xx);
                        plot(xx,yline,'-','Color',[0.6 0.6 0.6],'LineWidth',1.0);
                    catch, end
                end
                set(gca,'XColor',axisCol,'YColor',axisCol);


            else
                R=Pair(i,j).rho; P=Pair(i,j).p; rg=Pair(i,j).rho_g;
                txt=sprintf('Corr: %.3f%s', R, stars(P));
                if any(isfinite(rg))
                    t2=""; for g=1:numel(gNames), if isfinite(rg(g)), t2=t2+sprintf('\n%s: %.3f', gNames{g}, rg(g)); end, end
                    txt=txt+string(t2);
                end
                text(0.02,0.95,txt,'Units','normalized','VerticalAlignment','top','Color',fg,'FontWeight','bold');
                set(gca,'Visible','off');
            end
        else
            % --- 最右一列：该行变量按组箱线（兼容 R2024b） ---
            xi = X(:,i); mgAll = isfinite(xi) & isfinite(gLabels);
            if any(mgAll)
                hold on
                for gi = 1:numel(gNames)
                    mg = (gLabels == gi) & mgAll;
                    if any(mg)
                        bc = boxchart(gi * ones(sum(mg),1), xi(mg));  % 每组一个 box
                        bc.BoxFaceColor  = col(gi,:);                 % 组颜色
                        bc.BoxFaceAlpha  = 0.9;
                        bc.BoxEdgeColor  = [0.6 0.6 0.6];
                        bc.MarkerStyle   = 'none';                    % 不画散点
                        % 不设置 WhiskerLineStyle / WhiskerLineWidth（避免版本不兼容）
                    end
                end
                xlim([0.5, numel(gNames)+0.5]);
                set(gca,'XTick',1:numel(gNames),'XTickLabel',gNames,'XTickLabelRotation',45);
                set(gca,'YColor',axisCol);
            else
                set(gca,'Visible','off');
            end


        end
        if j>1 && j<=L, set(gca,'YTickLabel',[]); end
        if i<L  && j<=L, set(gca,'XTickLabel',[]); end
    end
end

str = sprintf('%s vs covariates (rank=%d, group=%s)', ...
      char(target), A.use_rank, char(A.group_mode));
hTitle = annotation(gcf,'textbox',[0.5 0.985 0 0], 'String',str, ...
    'FitBoxToText','on','Units','normalized', ...
    'HorizontalAlignment','center','VerticalAlignment','top', ...
    'LineStyle','none','Interpreter','none','FontWeight','bold');

if A.export_pdf
    set(fig,'Renderer','painters'); exportgraphics(fig, A.out_pdf, 'ContentType','vector');
    fprintf('[export] saved: %s\n', A.out_pdf);
end
if A.show_fig, figure(fig); drawnow; end   % 强制 FIG 弹窗

end

% ===== helpers（与上一个文件保持一致的白名单/映射/分组/配色） =====
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
        if isempty(u), gLabels=ones(height(T),1); gNames={'All'};
        else, gNames = cellstr(u); gLabels=zeros(height(T),1); for i=1:numel(u), gLabels(g==u(i))=i; end
        end
    case 'cluster'
        Z=X; for j=1:size(Z,2), Z(:,j)=(Z(:,j)-mean(Z(:,j),'omitnan'))/std(Z(:,j),'omitnan'); end
        m=all(isfinite(Z),2); Zc=Z(m,:); [idx,~]=kmeans(Zc, A.group_k,'Replicates',10,'Distance','sqeuclidean');
        gLabels=zeros(height(Z),1); gLabels(m)=idx; gLabels(~m)=mode(gLabels(gLabels>0));
        gNames=arrayfun(@(k)sprintf('Cluster%d',k),1:A.group_k,'UniformOutput',false);
    otherwise, error('group_mode 仅支持 auto/region/cluster');
end
end

function col = palette_colors(~, K)
base = [0 191 166; 58 134 255; 108 76 207]/255; if K<=size(base,1), col=base(1:K,:); else
t=linspace(1,size(base,1),K); col=interp1(1:size(base,1),base,t,'linear'); end
end

function s = stars(p)
if ~isfinite(p), s=''; return; end
if p<1e-3, s='***'; elseif p<1e-2, s='**'; elseif p<0.05, s='*'; else, s=''; end
end

function xw = winsorize_vec(x, plo, phi)
x = double(x); qlo=quantile(x,plo,'all'); qhi=quantile(x,phi,'all');
xw = x; xw(x<qlo)=qlo; xw(x>qhi)=qhi;
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
[~,ia]=unique([string({Map.col})' [Map.is_numeric]'],'rows','stable'); Map=Map(ia);
end
