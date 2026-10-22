function plot_china_choropleth_from_shp_with_inset(shpPath, mapCsvPath, valueVar, outPdf, nameField)
% 主幅（全国，无投影边框）+ 南海小窗（细矩形边框），绿色 30 级
% R2024b Mapping Toolbox

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

%% ---------- 绿色 30 级分级 ----------
Nbin = 30;
h = 130/360; s0=0.20; s1=0.90; v0=0.98; v1=0.55;
Cmap = hsv2rgb([h*ones(Nbin,1), linspace(s0,s1,Nbin)', linspace(v0,v1,Nbin)']);

edges = unique(quantile(v, linspace(0,1,Nbin+1)));
if numel(edges) < (Nbin+1), edges = linspace(min(v), max(v), Nbin+1); end
bin = nan(size(vals));
for i=1:numel(vals)
    if isfinite(vals(i))
        b = find(vals(i) >= edges(1:end-1) & vals(i) <= edges(2:end), 1, 'first');
        if isempty(b), b = Nbin; end
        bin(i) = b;
    end
end

%% ---------- 1) 主幅（去掉投影边框） ----------
fig = figure('Color','w'); set(fig,'Units','normalized');

ax_main = axes('Parent',fig,'Units','normalized',...
    'ActivePositionProperty','position','Position',[0.07 0.18 0.86 0.76]);
axes(ax_main);
axesm('lambert','MapLatLimit',[18 54],'MapLonLimit',[73 135], ...
      'Frame','off','Grid','off','MeridianLabel','off','ParallelLabel','off');  % 关键：Frame off
axis off; hold on;

for i=1:Nfeat
    if isempty(rings{i}), continue; end
    for k=1:numel(rings{i})
        la = rings{i}{k}.lat; lo = rings{i}{k}.lon;
        fc = isfinite(bin(i)) * Cmap(max(1,bin(i)),:) + ~isfinite(bin(i))*[0.85 0.85 0.85];
        patchm(la, lo, 0, 'FaceColor', fc, 'EdgeColor',[0.35 0.35 0.35], 'LineWidth',0.8);
    end
end
title(sprintf('%s — 30-level (green)', valueVar));
colormap(ax_main, Cmap);

cb = colorbar('peer',ax_main,'southoutside'); cb.Label.String = valueVar;
axpos = ax_main.Position; cb.Units='normalized';
cb.Position = [axpos(1)+0.05*axpos(3), max(0.06,axpos(2)-0.06), 0.90*axpos(3), 0.03];
try
    idx = round(linspace(1, Nbin, 6));
    cb.Ticks = (idx-0.5)/Nbin; cb.TickLabels = compose('%.1f', edges(idx));
catch
end

%% ---------- 2) 南海小窗（细矩形边框，投影边框关闭） ----------
latlimS = [3 40];  lonlimS = [105 123];  cutoffHN = 18.5;
ax_inset = axes('Parent',fig,'Units','normalized',...
    'ActivePositionProperty','position','Position',[0.73 0.23 0.22 0.28]);
axes(ax_inset);
axesm('lambert','MapLatLimit',latlimS,'MapLonLimit',lonlimS, ...
      'Frame','off','Grid','off','MeridianLabel','off','ParallelLabel','off');  % 关键：Frame off
axis off; hold on;

for i=1:Nfeat
    if isempty(rings{i}), continue; end
    for k=1:numel(rings{i})
        la = rings{i}{k}.lat; lo = rings{i}{k}.lon;
        if intersects_box(la,lo,latlimS,lonlimS) && max(la)<cutoffHN
            fc = isfinite(bin(i)) * Cmap(max(1,bin(i)),:) + ~isfinite(bin(i))*[0.90 0.90 0.90];
            patchm(la, lo, 0, 'FaceColor', fc, 'EdgeColor',[0.45 0.45 0.45], 'LineWidth',0.7);
        end
    end
end
setm(ax_inset,'MapLatLimit',latlimS,'MapLonLimit',lonlimS);
title('南海诸岛'); colormap(ax_inset, Cmap);

% 细矩形边框（annotation 按轴位置画矩形）
posInset = ax_inset.Position;   % [x y w h] 归一化到 figure
annotation(fig,'rectangle',posInset,'Color',[0.2 0.2 0.2],'LineWidth',0.8);

axes(ax_main);  % 回到主轴

%% ---------- 导出 ----------
if nargin>=4 && strlength(string(outPdf))>0
    exportgraphics(fig, outPdf, 'ContentType','vector', 'BackgroundColor','white');
end
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
