function plot_china_choropleth_from_shp_with_inset(shpPath, mapCsvPath, valueVar, outPdf, nameField, varargin)
% 主幅（全国，无投影边框）+ 南海小窗（细矩形边框，含上方大陆），绿色 30 级
% 并可选叠加“十段线/九段线”折线图层（经纬或投影坐标都可）。
%
% 必选参数：
%   shpPath     : 省级面 Shapefile（如 'T2024年初省级.shp'）
%   mapCsvPath  : CSV（含 {'province', valueVar}）
%   valueVar    : 'SRS' / 'LPI' 等
%   outPdf      : 导出 PDF（可空）
%   nameField   : 省名字段（推荐 'VAR_NAME'）
%
% Name-Value：
%   'UShapeShp'     : 断续线 Shapefile（如 '十段线.shp'，默认 '' 不叠加）
%   'DashColor'     : 断续线颜色（默认 [0.35 0.35 0.35]）
%   'DashWidth'     : 断续线线宽（默认 0.9）
%   'DashStyle'     : 断续线样式（默认 '--'）
%   'InsetLatLim'   : 小窗纬度范围（默认 [3 26.5]）
%   'InsetLonLim'   : 小窗经度范围（默认 [105 123]）
%   'HainanCutoff'  : 海南本岛最高纬度阈值（仅用于“仅画诸岛”情景，现已不再用，保留给兼容）

%% ---------- 解析可选参数 ----------
ip = inputParser;
addParameter(ip,'UShapeShp','',@(s)ischar(s)||isstring(s));
addParameter(ip,'DashColor',[0.35 0.35 0.35],@(v)isnumeric(v)&&numel(v)==3);
addParameter(ip,'DashWidth',0.9,@(x)isnumeric(x)&&isscalar(x));
addParameter(ip,'DashStyle','--',@(s)ischar(s)||isstring(s));
addParameter(ip,'InsetLatLim',[3 26.5],@(v)isnumeric(v)&&numel(v)==2);
addParameter(ip,'InsetLonLim',[105 123],@(v)isnumeric(v)&&numel(v)==2);
addParameter(ip,'HainanCutoff',18.5,@(x)isnumeric(x)&&isscalar(x)); %#ok<NVREPL> (兼容旧参数)
addParameter(ip,'ColorStyle','green',@(s)ischar(s)||isstring(s));
addParameter(ip,'NumLevels',30,@(x)isnumeric(x)&&isscalar(x)&&x>=3);

parse(ip,varargin{:});
P = ip.Results;

%% ---------- 读取省级面 ----------
S = shaperead(shpPath, 'UseGeoCoords', true);
hasLatLon = isfield(S,'Lat') && isfield(S,'Lon');
if ~hasLatLon
    S = shaperead(shpPath);
    assert(isfield(S,'X') && isfield(S,'Y'), 'Shapefile 坐标读取失败');
end

% 选用哪个属性做名称映射
allFields = fieldnames(S);
if nargin>=5 && strlength(string(nameField))>0
    useField = string(nameField);
else
    if any(strcmpi(allFields,'VAR_NAME')), useField="VAR_NAME";
    elseif any(strcmpi(allFields,'ENG_NAME')), useField="ENG_NAME";
    else, error('属性表缺少 VAR_NAME/ENG_NAME。可用字段：%s', strjoin(allFields,', '));
    end
end
fprintf('[map] using attribute: %s\n', useField);

%% ---------- 解析 rings 与省份键 ----------
Nfeat    = numel(S);
rings    = cell(Nfeat,1);
prov_key = strings(Nfeat,1);

for i = 1:Nfeat
    raw = string(S(i).(useField));
    if strcmpi(useField,'VAR_NAME')
        eng = ""; if any(strcmpi(allFields,'ENG_NAME')), eng = string(S(i).ENG_NAME); end
        prov_key(i) = var_name_to_key(raw, eng);
    elseif strcmpi(useField,'ENG_NAME')
        prov_key(i) = eng_name_to_key(raw);
    else
        prov_key(i) = cn_name_to_key(raw);
    end

    if hasLatLon, lat = double(S(i).Lat); lon = double(S(i).Lon);
    else,         lon = double(S(i).X);   lat = double(S(i).Y); end
    nanIdx = isnan(lat) | isnan(lon);
    cuts   = [0, find(nanIdx), numel(lat)+1];
    Ri = {};
    for k=1:numel(cuts)-1
        s = cuts(k)+1; e = cuts(k+1)-1;
        if e>=s
            la = lat(s:e); lo = lon(s:e);
            if numel(la)>=3
                if abs(la(1)-la(end))>1e-10 || abs(lo(1)-lo(end))>1e-10
                    la(end+1)=la(1); lo(end+1)=lo(1);
                end
                Ri{end+1} = struct('lat',la,'lon',lo); %#ok<AGROW>
            end
        end
    end
    rings{i} = Ri;
end

%% ---------- 读 CSV 并映射 ----------
T = readtable(mapCsvPath);
assert(ismember('province',T.Properties.VariableNames) && ismember(valueVar,T.Properties.VariableNames), ...
       'CSV 缺少 province 或 %s', valueVar);

T.province = canon_key(T.province);
prov_key   = canon_key(prov_key);

vals = nan(Nfeat,1);
for i=1:Nfeat
    if prov_key(i)~=""
        m = (T.province==prov_key(i));
        if any(m), vals(i) = T.(valueVar)(find(m,1)); end
    end
end
v = vals(isfinite(vals));  assert(~isempty(v),'没有可用的数值');


%% ---------- 可调色带 + 分级 ----------
Nbin = P.NumLevels;
% 优先用你的工厂函数；若不在 path，做一个最靠近的 HSV 回退
try
    Cmap = colormap_factory(P.ColorStyle, Nbin);
catch
    warning('colormap_factory 未找到，回退到 HSV 方案（style=%s, N=%d）。', string(P.ColorStyle), Nbin);
    % 简单回退：不同 style 映射到不同基色（保持与你的定义大致一致）
    switch lower(string(P.ColorStyle))
        case "orange", h = 28/360;  s0=0.25; s1=0.95; v0=1.00; v1=0.55;
        case "purple", h = 270/360; s0=0.20; s1=0.90; v0=0.98; v1=0.55;
        case "blue",   h = 210/360; s0=0.20; s1=0.95; v0=0.98; v1=0.55;
        case "red",    h =   5/360; s0=0.25; s1=0.95; v0=0.98; v1=0.55;
        case "yellow", h =  48/360; s0=0.25; s1=0.95; v0=1.00; v1=0.65;
        case "teal",   h = 170/360; s0=0.20; s1=0.85; v0=0.98; v1=0.60;
        otherwise      % 默认 green
            h = 130/360; s0=0.20; s1=0.90; v0=0.98; v1=0.55;
    end
    Cmap = hsv2rgb([h*ones(Nbin,1), linspace(s0,s1,Nbin)', linspace(v0,v1,Nbin)']);
end

% 分箱边界（分位数→去重→不足则线性）
edges = unique(quantile(v, linspace(0,1,Nbin+1)));
if numel(edges) < (Nbin+1), edges = linspace(min(v), max(v), Nbin+1); end

% 将每个省值映射到 1..Nbin 的箱
bin = nan(size(vals));
for i=1:numel(vals)
    if isfinite(vals(i))
        % 注意右端点归并：<= edges(2:end)
        b = find(vals(i) >= edges(1:end-1) & vals(i) <= edges(2:end), 1, 'first');
        if isempty(b), b = Nbin; end
        bin(i) = b;
    end
end


%% ---------- 1) 主幅（无投影边框） ----------
fig = figure('Color','w'); set(fig,'Units','normalized');

ax_main = axes('Parent',fig,'Units','normalized',...
    'ActivePositionProperty','position','Position',[0.07 0.18 0.86 0.76]);
axes(ax_main);
axesm('lambert','MapLatLimit',[18 54],'MapLonLimit',[73 135], ...
      'Frame','off','Grid','off','MeridianLabel','off','ParallelLabel','off');
axis off; hold on;

for i=1:Nfeat
    if isempty(rings{i}), continue; end
    for k=1:numel(rings{i})
        la = rings{i}{k}.lat; lo = rings{i}{k}.lon;
        fc = isfinite(bin(i)) * Cmap(max(1,bin(i)),:) + ~isfinite(bin(i))*[0.85 0.85 0.85];
        patchm(la, lo, 0, 'FaceColor', fc, 'EdgeColor',[0.35 0.35 0.35], 'LineWidth',0.8);
    end
end
title(sprintf('%s — 30-level', valueVar), fontsize=14);
colormap(ax_main, Cmap);

% 叠加“十段线”（若提供）
if strlength(string(P.UShapeShp))>0
    overlay_dash_layer(ax_main, P.UShapeShp, [18 54], [73 135], P.DashStyle, P.DashColor, P.DashWidth);
end

% 色标（仅跟主轴）
cb = colorbar('peer',ax_main,'southoutside'); cb.Label.String = valueVar;
axpos = ax_main.Position; cb.Units='normalized';
cb.Position = [axpos(1)+0.05*axpos(3), max(0.06,axpos(2)-0.06), 0.90*axpos(3), 0.03];
try
    idx = round(linspace(1, Nbin, 6));
    cb.Ticks = (idx-0.5)/Nbin; cb.TickLabels = compose('%.1f', edges(idx));
catch
end

%% ---------- 2) 南海小窗（含上方大陆；虚线边界；细矩形框） ----------
latlimS = P.InsetLatLim;  lonlimS = P.InsetLonLim;
ax_inset = axes('Parent',fig,'Units','normalized',...
    'ActivePositionProperty','position','Position',[0.66 0.18 0.15 0.21]);
axes(ax_inset);
axesm('lambert','MapLatLimit',latlimS,'MapLonLimit',lonlimS, ...
      'Frame','off','Grid','off','MeridianLabel','off','ParallelLabel','off');
axis off; hold on;

% 先填色（不描边），使虚线边界更清楚
for i=1:Nfeat
    if isempty(rings{i}), continue; end
    for k=1:numel(rings{i})
        la = rings{i}{k}.lat; lo = rings{i}{k}.lon;
        if intersects_box(la,lo,latlimS,lonlimS)
            fc = isfinite(bin(i)) * Cmap(max(1,bin(i)),:) + ~isfinite(bin(i))*[0.90 0.90 0.90];
            patchm(la, lo, 0, 'FaceColor', fc, 'EdgeColor','none');
        end
    end
end

% 多边形外轮廓 ≈ 海岸/国界（虚线）
edgeCol = [0.45 0.45 0.45];
for i=1:Nfeat
    if isempty(rings{i}), continue; end
    for k=1:numel(rings{i})
        la = rings{i}{k}.lat; lo = rings{i}{k}.lon;
        if intersects_box(la,lo,latlimS,lonlimS)
            plotm(la, lo, '--', 'Color', edgeCol, 'LineWidth', 0.8);
        end
    end
end

% 叠加“十段线”（若提供）
if strlength(string(P.UShapeShp))>0
    overlay_dash_layer(ax_inset, P.UShapeShp, latlimS, lonlimS, P.DashStyle, P.DashColor, P.DashWidth);
end

% 锁窗，标题 & 细矩形外框
setm(ax_inset,'MapLatLimit',latlimS,'MapLonLimit',lonlimS);
title('South China Sea'); colormap(ax_inset, Cmap);
posInset = ax_inset.Position;
annotation(fig,'rectangle',posInset,'Color',[0.2 0.2 0.2],'LineWidth',0.8);

axes(ax_main);  % 返回主轴

%% ---------- 导出 ----------
% if nargin>=4 && strlength(string(outPdf))>0
%     exportgraphics(fig, outPdf, 'ContentType','vector', 'BackgroundColor','white');
% end
end

% ================= helpers =================
function tf = intersects_box(lat, lon, latlim, lonlim)
minla = min(lat); maxla = max(lat);
minlo = min(lon); maxlo = max(lon);
tf = ~(maxla < latlim(1) || minla > latlim(2) || ...
       maxlo < lonlim(1) || minlo > lonlim(2));
end

function k = canon_key(s)
s = string(s);
s = lower(strrep(s," ",""));
s = remove_tones(s);
s = strrep(s,'ü','u');
s = regexprep(s,'huizu|zizhiqu|zhuangzu',''); % 保险去尾缀
k = s;
end

function key = var_name_to_key(varName, engName)
base = canon_key(strrep(string(varName)," ",""));
e = lower(string(engName));
if contains(base,"shanxi")
    if contains(e,"shaanxi"), key = "shaanxi"; return; end
    if contains(e,"shanxi"),  key = "shanxi";  return; end
end
key = base;
end

function key = eng_name_to_key(engName), key = canon_key(engName); end

function key = cn_name_to_key(cnName)
cnName = string(cnName);
mapCN = [
    "北京","beijing"; "天津","tianjin"; "河北","hebei"; "山西","shanxi"; "内蒙古","neimenggu";
    "辽宁","liaoning"; "吉林","jilin"; "黑龙江","heilongjiang"; "上海","shanghai"; "江苏","jiangsu";
    "浙江","zhejiang"; "安徽","anhui"; "福建","fujian"; "江西","jiangxi"; "山东","shandong";
    "河南","henan"; "湖北","hubei"; "湖南","hunan"; "广东","guangdong"; "广西","guangxi";
    "海南","hainan"; "重庆","chongqing"; "四川","sichuan"; "贵州","guizhou"; "云南","yunnan";
    "西藏","xizang"; "陕西","shaanxi"; "甘肃","gansu"; "青海","qinghai"; "宁夏","ningxia"; "新疆","xinjiang"
];
idx = find(mapCN(:,1)==cnName,1);
if isempty(idx), key = ""; else, key = mapCN(idx,2); end
end

function clean = remove_tones(pinyin)
clean = string(pinyin);
reps = {
    'ā','a'; 'á','a'; 'ǎ','a'; 'à','a';
    'ē','e'; 'é','e'; 'ě','e'; 'è','e';
    'ī','i'; 'í','i'; 'ǐ','i'; 'ì','i';
    'ō','o'; 'ó','o'; 'ǒ','o'; 'ò','o';
    'ū','u'; 'ú','u'; 'ǔ','u'; 'ù','u';
    'ǖ','ü'; 'ǘ','ü'; 'ǚ','ü'; 'ǜ','ü';
    'ń','n'; 'ň','n'; 'ǹ','n';
    'Ā','A'; 'Á','A'; 'Ǎ','A'; 'À','A';
    'Ē','E'; 'É','E'; 'Ě','E'; 'È','E';
    'Ī','I'; 'Í','I'; 'Ǐ','I'; 'Ì','I';
    'Ō','O'; 'Ó','O'; 'Ǒ','O'; 'Ò','O';
    'Ū','U'; 'Ú','U'; 'Ǔ','U'; 'Ù','U';
};
for r=1:size(reps,1)
    clean = strrep(clean, reps{r,1}, reps{r,2});
end
end

% ---- 读取线图层为经纬度段并绘图的封装 ----
function overlay_dash_layer(ax, lineShpPath, latlim, lonlim, style, color, width)
% 读取线要素 → 得到分段 Lat/Lon → 在给定 map-axes 上 plotm
Segs = read_line_latlon(lineShpPath);
if isempty(Segs), warning('未能从 %s 读取到线段；跳过断续线。', lineShpPath); return; end
axes(ax); hold on;
for s=1:numel(Segs)
    la = Segs{s}.lat; lo = Segs{s}.lon;
    if intersects_box(la,lo,latlim,lonlim)
        plotm(la, lo, style, 'Color', color, 'LineWidth', width);
    end
end
end

function Segs = read_line_latlon(lineShpPath)
% 返回 cell，元素为 struct('lat',[],'lon',[])
Segs = {};
try
    L = shaperead(lineShpPath,'UseGeoCoords',true);
    if isfield(L,'Lat') && isfield(L,'Lon')
        for i=1:numel(L)
            lat = double(L(i).Lat); lon = double(L(i).Lon);
            % 拆 NaN 分段
            cuts = [0, find(isnan(lat)|isnan(lon)), numel(lat)+1];
            for k=1:numel(cuts)-1
                s = cuts(k)+1; e = cuts(k+1)-1;
                if e>=s
                    Segs{end+1} = struct('lat',lat(s:e),'lon',lon(s:e)); %#ok<AGROW>
                end
            end
        end
        return;
    end
catch
    % 再尝试平面坐标读取
end

% 平面坐标 → 根据 PRJ 投影到经纬度
L = shaperead(lineShpPath);   % X/Y
assert(isfield(L,'X') && isfield(L,'Y'), '无法读取线要素 X/Y');
inf = shapeinfo(lineShpPath);
crs = [];
if isfield(inf,'CoordinateReferenceSystem')
    crs = inf.CoordinateReferenceSystem;
end
if ~isempty(crs) && ~(isa(crs,'map.geocrs'))
    % 投影坐标（projcrs），用 projinv 转回经纬
    for i=1:numel(L)
        x = double(L(i).X); y = double(L(i).Y);
        if all(isnan(x)|isnan(y)), continue; end
        [lat, lon] = projinv(crs, x, y);     % R2020b+ 可用
        cuts = [0, find(isnan(lat)|isnan(lon)), numel(lat)+1];
        for k=1:numel(cuts)-1
            s = cuts(k)+1; e = cuts(k+1)-1;
            if e>=s
                Segs{end+1} = struct('lat',lat(s:e),'lon',lon(s:e)); %#ok<AGROW>
            end
        end
    end
else
    % 没有 CRS 信息，就按经纬处理（很多断续线本来就是经纬）
    for i=1:numel(L)
        lon = double(L(i).X); lat = double(L(i).Y);
        cuts = [0, find(isnan(lat)|isnan(lon)), numel(lat)+1];
        for k=1:numel(cuts)-1
            s = cuts(k)+1; e = cuts(k+1)-1;
            if e>=s
                Segs{end+1} = struct('lat',lat(s:e),'lon',lon(s:e)); %#ok<AGROW>
            end
        end
    end
end
end
