function Tbl = screen_covariates(scores_root, features_mat, varargin)
% SCREEN_COVARIATES  计算目标指标与协变量的 Spearman 相关（含 p、FDR q）
% - 输出完整表（所有协变量），按 |rho| 从大到小排序
% - 落盘 CSV：scores/varscreen_<TARGET>.csv
% - 可选绘图：条形图（按 |rho| 降序），显著与否颜色区分
%
% 示例：
% Tbl = screen_covariates('scores','features_aug.mat', ...
%   'target','SRS','use_rank',true,'winsor',[0.01 0.99], ...
%   'alpha',0.05,'show_fig',true,'top_k_label',12);
%
% 返回表字段：
%   label        （白名单友好名）
%   column       （features 中实际列名）
%   n_eff        （有效样本量）
%   rho          （Spearman ρ）
%   p            （p 值）
%   q            （BH-FDR q 值）
%   signif       （q < alpha）
%   rank_abs     （|rho| 排名）

p = inputParser; p.FunctionName = mfilename;
addRequired(p,'scores_root',@(s)ischar(s)||isstring(s));
addRequired(p,'features_mat',@(s)ischar(s)||isstring(s));
addParameter(p,'target','SRS',@(s) any(strcmpi(string(s),{'RRI_VA','MEI_VA','SRS','LPI'})));
addParameter(p,'use_rank',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'winsor',[0.01 0.99],@(v)isnumeric(v)&&numel(v)==2&&v(1)>=0&&v(2)<=1&&v(1)<v(2));
addParameter(p,'alpha',0.05,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
addParameter(p,'show_fig',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'top_k_label',15,@(x)isnumeric(x)&&isscalar(x)&&x>=1); % 图上只标注前K名称，避免太挤
parse(p,scores_root,features_mat,varargin{:}); A=p.Results;

% ---------- 读取目标指标（省级） ----------
S1 = readtable(fullfile(A.scores_root,'map_ready_SRS.csv'),'VariableNamingRule','preserve');
S2 = readtable(fullfile(A.scores_root,'map_ready_SRS_LPI.csv'),'VariableNamingRule','preserve');
S1.province = string(S1.province); S2.province = string(S2.province);
K = innerjoin(S1(:,{'province','RRI_VA','MEI_VA','SRS'}), S2(:,{'province','LPI'}),'Keys','province');
K = sortrows(K,'province');
if A.use_rank
    for v = ["RRI_VA","MEI_VA","SRS","LPI"], K.(v)=tiedrank(double(K.(v))); end
end
target = upper(string(A.target));
if ~ismember(char(target), K.Properties.VariableNames)
    error('目标列 %s 不在指标表中。', char(target));
end

% ---------- 读取 features & 省级聚合 ----------
Feat = load_features_aug(A.features_mat);
if ~isstring(Feat.province), Feat.province = string(Feat.province); end
Feat = Feat(:, setdiff(Feat.Properties.VariableNames, {'wave_id'})); % 去掉逐波键
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

% ---------- 白名单（只评估你给的那 14 个；两项分类型跳过） ----------
covar_list = { ...
 'age65 share', ...
 'education share', ...
 'gdp_{pc}', ...
 'hospital beds per 1k', ...
 'proportion of students', ...
 'inbound overnight visitors', ...
 'disposable income_{pc}', ...
 'passenger volume', ...
 'public transit volume', ...
 'tertiary-sector share', ...
 'urban population density', ...
 'urban rate' };
% 这两项为分类型，此函数不做相关：'climate zone','region category'

% 映射到 features 列名（宽松匹配）
Map = map_covar_to_columns_local(Agg.Properties.VariableNames, covar_list);

% winsorize
AggW = Agg;
for k=1:height(Map)
    cn = Map.column{k};
    AggW.(cn) = winsorize_vec(double(AggW.(cn)), A.winsor(1), A.winsor(2));
end

% 合并到目标
T = innerjoin(K(:,["province",char(target)]), AggW, 'Keys','province');

% ---------- 逐变量相关 ----------
N = height(Map);
rho = nan(N,1); pval = nan(N,1); n_eff = zeros(N,1);
for i=1:N
    x = double(T.(Map.column{i}));
    y = double(T.(char(target)));
    m = isfinite(x) & isfinite(y);
    n_eff(i) = sum(m);
    if n_eff(i) >= 5
        [rho(i), pval(i)] = corr(x(m), y(m), 'Type','Spearman','Rows','complete');
    end
end
q = fdr_bh_local(pval);
[~,ord] = sort(abs(rho),'descend','MissingPlacement','last');

Tbl = table( string(Map.label), string(Map.column), n_eff, rho, pval, q, ...
    'VariableNames', {'label','column','n_eff','rho','p','q'});

% 先排序（|rho| 降序）
Tbl = Tbl(ord, :);

% 再计算名次（并列给平均名次）
Tbl.rank_abs = tiedrank(-abs(Tbl.rho));

Tbl.signif = Tbl.q < A.alpha;


% ---------- 落盘 ----------
outcsv = fullfile(A.scores_root, sprintf('varscreen_%s.csv', char(target)));
writetable(Tbl, outcsv);
fprintf('[screen_covariates] 全部协变量已保存：%s\n', outcsv);

% ---------- 可选：条形图 ----------
if A.show_fig
    fig = figure('Color','w','Position',[80 80 860 640]); movegui(fig,'center');
    ax = axes('Parent',fig); hold(ax,'on'); box(ax,'on');
    x = 1:height(Tbl);
    b = bar(ax, x, abs(Tbl.rho), 0.8, 'FaceColor','flat', 'EdgeColor',[0.2 0.2 0.2]);
    col = repmat([0.65 0.75 0.95],height(Tbl),1);
    col(Tbl.signif,:) = repmat([0.18 0.55 0.34], sum(Tbl.signif),1); % 显著=绿色
    b.CData = col;
    set(ax,'XLim',[0.5 height(Tbl)+0.5], 'XTick', x, ...
        'XTickLabel', string(Tbl.label), 'XTickLabelRotation', 45);
    ylabel(ax, '|Spearman \rho|'); title(ax, sprintf('Covariate screening — %s', char(target)));
    yline(ax,0,'k-');
    % 只标注前 K 个名称，避免太挤（已作为 XTickLabel，若想更简洁可只标注 top_k）
    if A.top_k_label < height(Tbl)
        keep = false(height(Tbl),1); keep(1:A.top_k_label)=true;
        xt = get(ax,'XTickLabel'); xt(~keep) = repmat("",sum(~keep),1);
        set(ax,'XTickLabel',xt);
    end
end
end

% ===== helpers =====

function Map = map_covar_to_columns_local(colNames, labels)
% 把给定的友好名映射到 features 列名（宽松，大小写/下划线均忽略）
C = string(colNames);
Ccanon = lower(regexprep(C,'[^a-z0-9]',''));
Map = table(strings(0,1), strings(0,1), 'VariableNames',{'label','column'});
for i=1:numel(labels)
    lab = string(labels{i});
    als = alias_of_label(lab);
    best = ""; bestScore = -inf;
    for j=1:numel(C)
        sc = -inf;
        for a = als
            tok = string(regexp(lower(regexprep(a,'[^a-z0-9]','')),'[a-z0-9]+','match'));
            if all(contains(Ccanon(j), tok))
                sc = max(sc, numel(tok) + 0.001*strlength(a));
            end
        end
        if sc > bestScore
            bestScore = sc; best = C(j);
        end
    end
    if bestScore > 0
        Map = [Map; {lab, best}]; %#ok<AGROW>
    else
        warning('未在 features 中匹配到：%s（跳过）', lab);
    end
end
% 去重
[~,ia] = unique(lower(Map.column),'stable'); Map = Map(ia,:);
end

function als = alias_of_label(lab)
switch lower(char(lab))
    case 'age65 share',                als = ["age65share","age65","elderlyshare","share65","pop65","age_65"];
    case 'education share',            als = ["educationshare","edu_share","tertiaryedu","college_share","higheredu"];
    case 'gdp_{pc}',                   als = ["gdppc","gdppercapita","gdp_pc","log_gdp_pc","gdp_percap"];
    case 'hospital beds per 1k',       als = ["hospitalbeds","bedsper1000","beds_per1k","beds_pc","bed_density"];
    case 'proportion of students',     als = ['proportion of students',"studentshare","proportionstudents","enrollment_share"];
    case 'inbound overnight visitors', als = ["overnightvisitors","inboundvisitors","tourist_night","visitor_night"];
    case 'disposable income_{pc}',     als = ["dispincpc","incomepercapita","disposable_income_pc","income_pc","log_income_pc"];
    case 'passenger volume',           als = ["passengervolume","passenger_flow","air_passenger","rail_passenger"];
    case 'public transit volume',      als = ["transitvolume","publictransit","log_transit_vol","metro_bus_rail"];
    case 'tertiary-sector share',      als = ["tertiaryshare","serviceshare","service_share","tertiary"];
    case 'urban population density',   als = ["urbanpopulationdensity","urban_density","pop_density_urban"];
    case 'urban rate',                 als = ["urbanrate","urbanization","urban_ratio","urban_share"];
    otherwise,                         als = string(lab);
end
end

function xw = winsorize_vec(x, plo, phi)
x = double(x); qlo = quantile(x,plo,'all'); qhi = quantile(x,phi,'all');
xw = x; xw(x<qlo)=qlo; xw(x>qhi)=qhi;
end

function q = fdr_bh_local(p)
q = nan(size(p)); ix=find(isfinite(p)); if isempty(ix), return; end
[sp,ord]=sort(p(ix)); m=numel(sp);
qv = sp .* m ./ (1:m)'; qv = flipud(cummin(flipud(qv)));
q(ix(ord)) = qv;
end
