function plot_waves_overlay_with_insets(varargin)
% 单轴叠加（不归一化）+ 指定三大波 inset：
%   A: 2020–2021，嵌 Hubei 第1波
%   B: 2022       嵌 Shanghai 第2波、Jilin 第3波
%
% 目录：<root>/<Province>/ResultsXX/daily_cases.xlsx
% 元数据：features_aug.mat 内表 Features 必含 {province, wave_id('Hubei-01'等), start_date}
%
% ===== 参数 =====
ip = inputParser;
addParameter(ip,'rootDir','.',@(s)ischar(s)||isstring(s));
addParameter(ip,'featuresMat','features_aug.mat',@(s)ischar(s)||isstring(s));
addParameter(ip,'outA','Fig_Overlay_2020_2021',@(s)ischar(s)||isstring(s));
addParameter(ip,'outB','Fig_Overlay_2022',@(s)ischar(s)||isstring(s));
% 风格 & 字号
addParameter(ip,'yScale','linear',@(s)any(strcmpi(s,{'linear','log'})));
addParameter(ip,'smoothDays',7,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(ip,'FontName','Helvetica',@(s)ischar(s)||isstring(s));
addParameter(ip,'FS_Title',10,@(x)isnumeric(x)&&isscalar(x));
addParameter(ip,'FS_Axis',8,@(x)isnumeric(x)&&isscalar(x));
addParameter(ip,'FS_Legend',8,@(x)isnumeric(x)&&isscalar(x));
addParameter(ip,'LineWidth',2,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'LineLighten',0.12,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=0.9);
% 阶段背景带
addParameter(ip,'ShowPhases',true,@islogical);
% 嵌入或分开
addParameter(ip,'EmbedInsets',true,@islogical);

addParameter(ip,'ExcludeWaves', {'Hubei',1; 'Shanghai',2; 'Jilin',3});  % {province, wave_no}

addParameter(ip,'LegendRows',2,@(x)isnumeric(x)&&isscalar(x)&&x>=1);   % 图例行数：2 或 3
addParameter(ip,'LegendToken',16,@(x)isnumeric(x)&&isscalar(x)&&x>0);  % 色标小方块尺寸(px)


parse(ip,varargin{:});
P = ip.Results;

%% --- 时段与切分 ---
rngA = [datetime(2020,1,1) datetime(2021,12,31)];
rngB = [datetime(2022,1,1) datetime(2022,10,31)];
cut1 = datetime(2021,6,1);  % Original/Alpha → Delta
cut2 = datetime(2022,1,1);  % Delta → Omicron

%% --- 读取 Features & 波数据缓存（一次读，复用两图） ---
S = load(char(P.featuresMat));
if isfield(S,'Features'), T=S.Features; else
    fns=fieldnames(S); T=[]; for i=1:numel(fns), if istable(S.(fns{i})), T=S.(fns{i}); break; end, end
end
assert(~isempty(T),'未在 %s 找到表。', P.featuresMat);

vnames = lower(string(T.Properties.VariableNames));
province_col = T.Properties.VariableNames{find(vnames=="province",1)};
waveid_col   = T.Properties.VariableNames{find(vnames=="wave_id"|vnames=="waveid",1)};
start_col    = T.Properties.VariableNames{find(vnames=="start_date"|vnames=="startdate",1)};
rid = find(vnames=="region category"|vnames=="region_category"|vnames=="region",1);
assert(~isempty(province_col)&&~isempty(waveid_col)&&~isempty(start_col), 'Features 缺列。');
region_col = ''; if ~isempty(rid), region_col=T.Properties.VariableNames{rid}; end

if ~isdatetime(T.(start_col))
    try, T.(start_col)=datetime(T.(start_col));
    catch, T.(start_col)=datetime(T.(start_col),'InputFormat','yyyy-MM-dd');
    end
end
T.(province_col)=string(T.(province_col));
T.(waveid_col)=string(T.(waveid_col));
if ~isempty(region_col), T.(region_col)=string(T.(region_col)); end

% 排序（区域块内字母序），保证上下两图颜色映射一致
if isempty(region_col)
    provOrder = sort(unique(T.(province_col),'stable'));
else
    regs = unique(T.(region_col),'stable');
    provOrder = strings(0,1);
    for i=1:numel(regs)
        p = sort(unique(T.(province_col)(T.(region_col)==regs(i)),'stable'));
        provOrder=[provOrder; p]; %#ok<AGROW>
    end
end
Nprov = numel(provOrder);
[~, w_no] = parse_wave_id(T.(waveid_col)); T.wave_no = w_no;

% 读波数据（原始量级）
cache = table('Size',[height(T) 3], ...
    'VariableTypes',{'string','double','cell'}, ...
    'VariableNames',{'province','wave_no','series'});
fprintf('[scan] reading daily_cases...\n');
for i=1:height(T)
    prov = char(T.(province_col)(i)); wn = T.wave_no(i);
    if isnan(wn), continue; end
    xls = fullfile(char(P.rootDir), prov, sprintf('Results%02d', wn), 'daily_cases.xlsx');
    if ~isfile(xls), warning('缺少 %s（跳过）', xls); continue; end
    [dt,y] = read_cases_xlsx(xls);
    if isempty(dt)||isempty(y), warning('空数据：%s', xls); continue; end
    if P.smoothDays>1, y = movmean(y,[floor((P.smoothDays-1)/2) ceil((P.smoothDays-1)/2)],'omitnan'); end
    cache.province(i)=string(prov);
    cache.wave_no(i)=wn;
    cache.series{i}=table(dt(:),y(:),'VariableNames',{'date','cases'});
end

% 颜色：30 色自然风（固定映射）
C30 = lighten_colors(palette_nature30(max(Nprov,30)), P.LineLighten);
provColor = containers.Map('KeyType','char','ValueType','any');
for i=1:Nprov, provColor(char(provOrder(i))) = C30(i,:); end

%% --- 画图 A（2020–2021） ---
draw_main_panel(cache, provOrder, provColor, rngA(1), rngA(2), P, P.outA, P.ShowPhases, [cut1 cut2], P.ExcludeWaves,'2020');
% 在 A 中嵌 Hubei 第1波
if P.EmbedInsets
    draw_one_inset(gcf, gca, cache, provColor, 'Hubei', 1, rngA, P, [0.14 0.62 0.25 0.26]); % 左上合适处
else
    draw_standalone_inset(cache, provColor, 'Hubei', 1, rngA, P, 'Inset_Hubei_01');
end

%% --- 画图 B（2022） ---
draw_main_panel(cache, provOrder, provColor, rngB(1), rngB(2), P, P.outB, P.ShowPhases, [cut1 cut2], P.ExcludeWaves,'2022');
% 在 B 中嵌 Shanghai-02, Jilin-03
if P.EmbedInsets
    draw_one_inset(gcf, gca, cache, provColor, 'Shanghai', 2, rngB, P, [0.10 0.62 0.26 0.30]); % 左上
    draw_one_inset(gcf, gca, cache, provColor, 'Jilin',    3, rngB, P, [0.38 0.62 0.26 0.30]); % 中上
else
    draw_standalone_inset(cache, provColor, 'Shanghai', 2, rngB, P, 'Inset_Shanghai_02');
    draw_standalone_inset(cache, provColor, 'Jilin',    3, rngB, P, 'Inset_Jilin_03');
end
end

%% ========= 子程序 =========
function draw_main_panel(cache, provOrder, provColor, xmin, xmax, P, outName, showPhases, cuts, excludeWaves,nianfen)
fig = figure('Color','w','Name',char(outName),'Units','pixels','Position',[80 80 1400 900]);
ax = axes('Parent',fig); hold(ax,'on'); box(ax,'off');
set(ax,'FontName',P.FontName,'FontSize',P.FS_Axis,'LineWidth',0.9);
% 时间轴
if nianfen == '2020'
    setup_time_axis(ax, xmin, xmax, P.FS_Axis,2);
elseif nianfen == '2022'
    setup_time_axis(ax, xmin, xmax, P.FS_Axis,1);
end
% 背景带
if showPhases
    draw_phase_bands(ax, xmin, xmax, cuts);
end
% 叠加所有省所有波
for ip=1:numel(provOrder)
    prov = char(provOrder(ip)); col = provColor(prov);
    rows = find(cache.province==string(prov));
    for r=rows'
        if is_excluded(prov, cache.wave_no(r), excludeWaves), continue; end
        S = cache.series{r}; if isempty(S), continue; end
        [d1,y1] = mask_range(S.date, S.cases, xmin, xmax);
        if numel(d1)<2 || all(~isfinite(y1)), continue; end
        if strcmpi(P.yScale,'log'), yplot = log10(max(0,y1)+1); else, yplot=y1; end
        plot(ax, datenum(d1), yplot, '-', 'Color', col, 'LineWidth', P.LineWidth, 'HandleVisibility','off');
    end
end
% 轴标签

% ---- Y 轴标签 ----
if strcmpi(P.yScale,'log')
    yl = "Daily cases (log_{10}(y+1))";
else
    yl = "Daily cases";
end
ylabel(ax, yl, 'FontSize', P.FS_Axis);
% xlabel(ax,'Date');
title(ax, 'COVID-19 in China (2022)', 'FontWeight','bold','FontSize',P.FS_Title);

% xlabel(ax,'Date');

title(ax, 'COVID-19 in China (2020–2021)', 'FontWeight','bold','FontSize',P.FS_Title);


% ------- Bottom legend: 30 省横向长条（2~3 行） -------
% 1) 先为每省做一根“代表线”，仅用于图例
hLeg = gobjects(numel(provOrder),1);
for i = 1:numel(provOrder)
    c = provColor(char(provOrder(i)));
    hLeg(i) = plot(ax, NaN, NaN, '-', 'Color', c, 'LineWidth', P.LineWidth, ...
                   'HitTest','off');
end

% 2) 计算列数：按行数平均分布
nItems = numel(provOrder);
nCols  = ceil(nItems / max(1, P.LegendRows));

% 3) 生成底部横向图例
lgd = legend(ax, hLeg, cellstr(provOrder), ...
    'Location','southoutside', 'Orientation','horizontal', ...
    'NumColumns', nCols, 'Box','off', ...
    'Interpreter','none', 'FontName', P.FontName, 'FontSize', P.FS_Legend);

% 4) 图例外观微调
lgd.ItemTokenSize = [P.LegendToken  P.LegendToken];   % 色块大小
lgd.TextColor     = [0.1 0.1 0.1];

% 5) 给底部图例留一点空间（把主轴上移一点）
pos = ax.Position;
ax.Position = [pos(1), pos(2)+0.04, pos(3), pos(4)-0.06];  % 视需要再微调


% 右侧省名色签（两列）
% draw_sidepanel_prov_legend(fig, ax, provOrder, provColor, P.FS_Legend, P.FontName, 2);
% 导出
try
    savefig(fig, [char(outName) '.fig']);
    exportgraphics(fig,[char(outName) '.pdf'],'ContentType','vector','BackgroundColor','white');
catch, end
end

function draw_one_inset(fig, axMain, cache, provColor, provName, waveNo, trange, P, normPos)
% 在主图里嵌入一个 axes 做 inset。normPos=[x y w h] (normalized)
axes(axMain);
ax = axes('Parent',fig,'Position',normPos,'Box','on','LineWidth',0.6,'FontName',P.FontName);
hold(ax,'on');
set(ax,'FontSize',P.FS_Axis-0.5);
col = provColor(provName);
% 找该省该波的 series（若多条匹配则取 waveNo）
rowsAll = find(cache.province==string(provName));
row = NaN;
for r = rowsAll', if cache.wave_no(r)==waveNo, row=r; break; end, end
if isnan(row), warning('未找到 %s-%02d，小窗略过。', provName, waveNo); delete(ax); return; end
S = cache.series{row};
[d1,y1] = mask_range(S.date, S.cases, trange(1), trange(2));
if isempty(d1), warning('%s-%02d 在该时段无数据。', provName, waveNo); delete(ax); return; end
if strcmpi(P.yScale,'log'), yplot = log10(max(0,y1)+1); else, yplot = y1; end
plot(ax, datenum(d1), yplot, '-', 'Color', col, 'LineWidth', max(1.3,P.LineWidth+0.4));
% 美化：极简轴、少刻度
setup_time_axis(ax, min(d1), max(d1), P.FS_Axis-0.5);
grid(ax,'off'); ax.YGrid='on'; ax.GridAlpha=0.08;
title(ax, sprintf('%s (%s)', provName, choose_label(trange)), 'FontWeight','normal','FontSize',P.FS_Axis+0.2,'Color',[0.1 0.1 0.1]);
end

function draw_standalone_inset(cache, provColor, provName, waveNo, trange, P, outName)
fig = figure('Color','w','Name',char(outName),'Units','pixels','Position',[200 200 600 360]);
ax = axes('Parent',fig); hold(ax,'on'); box(ax,'on'); set(ax,'FontName',P.FontName,'FontSize',P.FS_Axis);
col = provColor(provName);
rowsAll = find(cache.province==string(provName));
row = NaN; for r=rowsAll', if cache.wave_no(r)==waveNo, row=r; break; end, end
if isnan(row), warning('未找到 %s-%02d。', provName, waveNo); return; end
S = cache.series{row}; [d1,y1] = mask_range(S.date, S.cases, trange(1), trange(2));
if isempty(d1), warning('%s-%02d 在该时段无数据。', provName, waveNo); return; end
if strcmpi(P.yScale,'log'), yplot=log10(max(0,y1)+1); else, yplot=y1; end
plot(ax, datenum(d1), yplot, '-', 'Color', col, 'LineWidth', max(1.3,P.LineWidth+0.4));
setup_time_axis(ax, min(d1), max(d1), P.FS_Axis);
title(ax, sprintf('%s (%s)', provName, choose_label(trange)), 'FontWeight','normal','FontSize',P.FS_Axis+0.5);
try
    savefig(fig, [char(outName) '.fig']); exportgraphics(fig,[char(outName) '.pdf'],'ContentType','vector','BackgroundColor','white');
catch, end
end

function label = choose_label(trange)
y = year(trange(1));
if y==2020||y==2021, label='Original/Delta'; else, label='Omicron'; end
end

%% --- 共用工具 ---
function [d,y]=read_cases_xlsx(xlsfile)
T=readtable(xlsfile,'VariableNamingRule','preserve');
if size(T,1)<2||size(T,2)<2, d=datetime.empty; y=[]; return; end
d=T{2:end,1}; if ~isdatetime(d), try, d=datetime(d); catch, d=datetime(d,'InputFormat','yyyy-M-d'); end, end
y=double(T{2:end,2});
end

function [prov,no]=parse_wave_id(wave_id)
wid=string(wave_id); prov=strings(numel(wid),1); no=zeros(numel(wid),1);
for i=1:numel(wid)
    s=char(wid(i)); k=strfind(s,'-');
    if isempty(k), prov(i)=string(s); no(i)=NaN;
    else
        prov(i)=string(s(1:k(end)-1));
        numstr=regexprep(s(k(end)+1:end),'[^\d]',''); v=str2double(numstr); if isnan(v), v=NaN; end
        no(i)=v;
    end
end
end

function setup_time_axis(ax, xmin, xmax, FS, tickMonths)
xlim(ax, datenum([xmin xmax]));
ax.XLimMode='manual';
if nargin<5 || isempty(tickMonths), tickMonths = 2; end   % 默认隔月
xt = datenum(dateshift(xmin,'start','month'):calmonths(tickMonths):dateshift(xmax,'start','month'));
set(ax,'XTick',xt);
datetick(ax,'x','yyyy-mmm','keepticks','keeplimits');
ax.XTickMode='manual';
grid(ax,'on'); ax.GridAlpha=0.10; ax.YGrid='on'; ax.MinorGridAlpha=0.08; ax.Layer='top';
ax.FontSize=FS; ax.YAxis.Color=[0.1 0.1 0.1];
end

function [d2,y2]=mask_range(d,y,xmin,xmax)
m=(d>=xmin)&(d<=xmax)&isfinite(y);
d2=d(m); y2=y(m);
end

function draw_phase_bands(ax, xmin, xmax, cuts)
yl=ylim(ax); ylo=yl(1); yhi=yl(2);
XA=datenum([xmin, min(xmax, cuts(1)-days(1))]);
XD=datenum([cuts(1), cuts(2)-days(1)]);
XO=datenum([cuts(2), xmax]);
colA=[0.94 0.94 0.94]; colD=[0.92 0.92 0.96]; colO=[0.95 0.92 0.97];
if XA(2)>XA(1), patch('XData',[XA(1) XA(2) XA(2) XA(1)],'YData',[ylo ylo yhi yhi],'FaceColor',colA,'EdgeColor','none','Parent',ax,'HitTest','off'); end
if XD(2)>XD(1), patch('XData',[XD(1) XD(2) XD(2) XD(1)],'YData',[ylo ylo yhi yhi],'FaceColor',colD,'EdgeColor','none','Parent',ax,'HitTest','off'); end
if XO(2)>XO(1), patch('XData',[XO(1) XO(2) XO(2) XO(1)],'YData',[ylo ylo yhi yhi],'FaceColor',colO,'EdgeColor','none','Parent',ax,'HitTest','off'); end
uistack(findobj(ax,'Type','line'),'top');
end

function draw_sidepanel_prov_legend(fig, ax, provOrder, provColor, FS, fontName, nCols)
pos = ax.Position; margin=0.02; panelWidth=0.20;
ax.Position = [pos(1) pos(2) (1-panelWidth-2*margin - pos(1)) pos(4)];
panelLeft = ax.Position(1)+ax.Position(3)+margin;
uip = uipanel('Parent',fig,'Units','normalized','Position',[panelLeft ax.Position(2) panelWidth ax.Position(4)], ...
              'BorderType','none','BackgroundColor','w');
uicontrol('Style','text','Parent',uip,'Units','normalized','Position',[0 0.94 1 0.05], ...
          'String','Provinces','FontName',fontName,'FontSize',FS,'FontWeight','bold', ...
          'BackgroundColor','w','HorizontalAlignment','left');
N=numel(provOrder); nRows=ceil(N/nCols); rowH=0.90/nRows; colW=1.0/nCols;
for i=1:N
    col = provColor(char(provOrder(i)));
    c = mod(i-1,nCols); r=floor((i-1)/nCols);
    x=c*colW; y=0.90-(r+1)*rowH;
    uipanel('Parent',uip,'Units','normalized','Position',[x+0.02 y+rowH*0.28 0.035 rowH*0.44], ...
            'BackgroundColor',col,'BorderType','line','HighlightColor',[0.85 0.85 0.85],'ShadowColor',[0.85 0.85 0.85]);
    uicontrol('Style','text','Parent',uip,'Units','normalized','Position',[x+0.06 y colW-0.07 rowH], ...
              'String',char(provOrder(i)),'FontName',fontName,'FontSize',FS, ...
              'BackgroundColor','w','HorizontalAlignment','left');
end
end

function C = palette_nature30(N)
% 自然期刊风 30 色（避免刺眼黄，控制明度/饱和度）
if nargin<1, N=30; end
C=zeros(N,3); phi=137.508; h0=15;
for i=1:N
    h = mod(h0+(i-1)*phi,360);
    s = 0.55 + 0.10*sin((i-1)*pi/6);
    v = 0.78 + 0.08*cos((i-1)*pi/8);
    C(i,:) = hsv2rgb([h/360 s v]);
end
C=min(max(C,0),1);
end

function C2 = lighten_colors(C, amt)
if amt<=0, C2=C; else, C2 = C*(1-amt) + amt; end
C2=min(max(C2,0),1);
end



function tf = is_excluded(prov, wn, excl)
tf = false;
if isempty(excl), return; end
if istable(excl)
    provs = string(excl{:,1});
    waves = double(excl{:,2});
elseif iscell(excl)
    provs = string(excl(:,1));
    waves = cellfun(@double, excl(:,2));
else
    return;
end
tf = any(strcmpi(provs, string(prov)) & (waves(:) == wn));
end
