function plot_waves_210_single(varargin)
% PLOT_WAVES_210_SINGLE
% China’s 210 outbreak waves across 30 provinces (2020–2022) — 单轴岭线图
% 目录: <root>/<Province>/ResultsXX/daily_cases.xlsx
% features_aug.mat 的 Features 表需含: province, wave_id('Anhui-06'), start_date
% 可选: variant, region / region_category / region category

%% -------- 参数 --------
ip = inputParser;
addParameter(ip,'rootDir','.',@(s)ischar(s)||isstring(s));
addParameter(ip,'featuresMat','features_aug.mat',@(s)ischar(s)||isstring(s));
addParameter(ip,'outName','Fig_Waves210_single',@(s)ischar(s)||isstring(s));
addParameter(ip,'smoothDays',7,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(ip,'leftMarginDays',100,@(x)isnumeric(x)&&isscalar(x)&&x>0);   % 左侧留白（天数）
addParameter(ip,'labelShiftDays',110,@(x)isnumeric(x)&&isscalar(x));        % 文字再往左移多少天

parse(ip,varargin{:});
P = ip.Results;
rootDir=char(P.rootDir); featuresMat=char(P.featuresMat); outName=char(P.outName);
kSmooth=P.smoothDays;

leftMarginDays = P.leftMarginDays;
labelShiftDays = P.labelShiftDays;


% 时间范围与关键分界
xmin = datetime(2020,1,1); xmax = datetime(2022,12,31);
cut1 = datetime(2021,6,1);  % Original/Alpha → Delta
cut2 = datetime(2022,1,1);  % Delta → Omicron
leftMarginDays = 25;        % 给左侧省名留白

% 字体与样式
FS_title=16; FS_axis=13; FS_label=14;
alphaFill = 0.48;

%% -------- 读 Features --------
S = load(featuresMat);
T = []; if isfield(S,'Features'), T=S.Features; else
    fns=fieldnames(S); for i=1:numel(fns), if istable(S.(fns{i})), T=S.(fns{i}); break; end, end
end
assert(~isempty(T),'未在 %s 中找到表。', featuresMat);

vnames = lower(string(T.Properties.VariableNames));
province_col = T.Properties.VariableNames{find(vnames=="province",1)};
waveid_col   = T.Properties.VariableNames{find(vnames=="wave_id"|vnames=="waveid",1)};
start_col    = T.Properties.VariableNames{find(vnames=="start_date"|vnames=="startdate",1)};
rid = find(vnames=="region category" | vnames=="region_category" | vnames=="region",1);
assert(~isempty(province_col)&&~isempty(waveid_col)&&~isempty(start_col), ...
    'Features 缺少 province / wave_id / start_date');
region_col = ''; if ~isempty(rid), region_col=T.Properties.VariableNames{rid}; end

if ~isdatetime(T.(start_col))
    try, T.(start_col)=datetime(T.(start_col)); ...
    catch, T.(start_col)=datetime(T.(start_col),'InputFormat','yyyy-MM-dd'); end
end
T.(province_col)=string(T.(province_col));
T.(waveid_col)=string(T.(waveid_col));
if ~isempty(region_col), T.(region_col)=string(T.(region_col)); end

% 行排序：按 7 大区分块（有则用），块内按省名首字母序
if isempty(region_col)
    provOrder = sort(unique(T.(province_col),'stable'));
    blocks={provOrder}; blockNames={''};
else
    regs = unique(T.(region_col),'stable');
    provOrder=strings(0,1); blocks=cell(numel(regs),1); blockNames=cell(numel(regs),1);
    for i=1:numel(regs)
        p = sort(unique(T.(province_col)(T.(region_col)==regs(i)),'stable'));
        blocks{i}=p; blockNames{i}=char(regs(i)); provOrder=[provOrder; p]; %#ok<AGROW>
    end
end

% === 颜色：7 区域从上到下的“日照金山” ===
numRegs = max(1, numel(blockNames));              % 期望 7
regColors = ridge_row_cmap_local(numRegs);        % 顶金 → 底蓝

% === 为每个省标记其区域索引（用于着色） ===
regIdxPerProv = ones(numel(provOrder),1);
if ~isempty(region_col)
    regs = unique(T.(region_col),'stable');       % 与 blocks/blockNames 同顺序
    for i = 1:numel(provOrder)
        p = provOrder(i);
        ridx = find(regs == T.(region_col)(find(T.(province_col)==p,1,'first')), 1);
        if ~isempty(ridx), regIdxPerProv(i) = ridx; end
    end
end


Nprov = numel(provOrder); assert(Nprov>0,'没有省份可绘制');

% wave_id → 号
[~, w_no] = parse_wave_id(T.(waveid_col)); T.wave_no = w_no;

%% -------- 预读 & 省内峰值 --------
provMax = containers.Map('KeyType','char','ValueType','double');
cache = table('Size',[height(T) 5], ...
    'VariableTypes',{'string','double','datetime','double','cell'}, ...
    'VariableNames',{'province','wave_no','t0','peak','series'});

fprintf('[scan] reading waves...\n');
for i=1:height(T)
    prov = char(T.(province_col)(i)); wn = T.wave_no(i);
    if isnan(wn), continue; end
    xls = fullfile(rootDir, prov, sprintf('Results%02d', wn), 'daily_cases.xlsx');
    if ~isfile(xls), warning('缺少文件：%s（跳过）', xls); continue; end
    try
        [dt,y] = read_cases_xlsx(xls);
        if isempty(dt)||isempty(y), warning('空数据：%s', xls); continue; end
        if kSmooth>1, y=movmean(y,[floor((kSmooth-1)/2) ceil((kSmooth-1)/2)],'omitnan'); end
        cache.province(i)=string(prov); cache.wave_no(i)=wn; cache.t0(i)=dt(1);
        cache.series{i}=table(dt(:),y(:),'VariableNames',{'date','cases'});
        cache.peak(i)=max(y,[],'omitnan');
        if ~isKey(provMax,prov), provMax(prov)=0; end
        provMax(prov) = max(provMax(prov), cache.peak(i));
    catch ME
        warning('读取失败 %s：%s', xls, ME.message);
    end
end
kP=provMax.keys; for k=1:numel(kP), if provMax(kP{k})<=0, provMax(kP{k})=1; end, end

%% -------- 画布与单轴 --------
fig = figure('Color','w','Name','Waves210 SINGLE','Units','pixels','Position',[80 80 1400 900]);
ax = axes('Parent',fig); hold(ax,'on'); box(ax,'off');

% 行距与行配色
rowGap=0.1; rowH=1.0; rowStep=rowH+rowGap;
rowC = ridge_row_cmap_local(Nprov);
rowLine = min(rowC*0.85 + 0.15, 1); rowFill = rowC;

% 轴设置（统一 datenum）

setup_time_axis(ax, xmin - days(leftMarginDays), xmax, FS_axis);



% 左侧省名与区域标签

draw_left_labels(ax, provOrder, blocks, blockNames, ...
                 rowStep, FS_label, leftMarginDays, labelShiftDays, xmin, ...
                 regIdxPerProv, regColors);




% 分界竖线（同一轴）
vline(ax, cut1, 0.92, [0.55 0.55 0.55], ':', 0.9);
vline(ax, cut2, 0.92, [0.55 0.55 0.55], ':', 0.9);

%% -------- 绘制所有省的岭线 --------
fprintf('[plot] drawing ridgelines (single axis)...\n');
for ip=1:Nprov
    prov = char(provOrder(ip));
    y0 = (Nprov - ip) * rowStep;
    scale = 1 / provMax(prov);

    rows = find(cache.province==string(prov));
    if isempty(rows), continue; end
    peaks = cache.peak(rows);
    v = peaks(isfinite(peaks));
    if isempty(v), thr50=inf; thr75=inf; thr90=inf;
    else, q=quantile(v,[0.5 0.75 0.9]); thr50=q(1); thr75=q(2); thr90=q(3);
    end

    for r=rows'
        Srow=cache.series{r}; if isempty(Srow), continue; end
        d=Srow.date; y=Srow.cases; if all(~isfinite(y)), continue; end
        yN = y*scale;
        [d1,y1] = mask_range(d,yN,xmin,xmax);
        if numel(d1)>=2 && any(isfinite(y1))
            lw = choose_lw(cache.peak(r),thr50,thr75,thr90);
            draw_one(ax, d1, y1, y0, rowLine(ip,:), rowFill(ip,:), alphaFill, lw);
        end
    end

    % % 极端峰值提示（>省内95%）
    % pk95 = safe_quantile(peaks,0.95);
    % [pkVal, idxMax] = max(peaks,[],'omitnan');
    % if pkVal>pk95 && ~isempty(idxMax)
    %     rr = rows(idxMax); SS = cache.series{rr};
    %     [~,imax]=max(SS.cases);
    %     dpk=SS.date(imax); ypk=(SS.cases(imax)*scale)+y0;
    %     if dpk>=xmin && dpk<=xmax
    %         text(ax, datenum(dpk), ypk, sprintf('  %.0f', pkVal), ...
    %             'VerticalAlignment','bottom','FontSize',FS_axis,'Color',[0.25 0.25 0.25]);
    %     end
    % end
end

% 外观
ylim(ax,[-0.20, Nprov*rowStep-0.05]);
set(ax,'YColor','none','TickDir','out','LineWidth',0.9,'FontName','Helvetica','FontSize',FS_axis);
xlabel(ax,'Date',FontSize=13); title(ax,"China’s 210 outbreak waves across 30 provinces (2020–2022)", ...
    'FontWeight','bold','FontSize',FS_title);

% 归一化注记
txt = "Height scaled to each province's 2020–2022 maximum (0–1); curves are not cross-province comparable.";
text(ax, datenum(xmin + days(6-leftMarginDays)), -0.12, txt, 'FontSize',FS_axis, ...
    'Color',[0.35 0.35 0.35],'VerticalAlignment','top');

datetick(ax,'x','yyyy-mmm','keepticks','keeplimits');  % ← 追加这一行

%% -------- 导出 FIG/PDF --------
figFile=[outName '.fig']; pdfFile=[outName '.pdf'];
try, savefig(fig,figFile); exportgraphics(fig,pdfFile,'ContentType','vector','BackgroundColor','white');
    fprintf('[out] saved: %s, %s\n', figFile, pdfFile);
catch ME, warning(sprintf('导出失败：%s', ME.message));
end
end

%% ===== 工具函数 =====
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
function setup_time_axis(ax,xmin,xmax,FS)
xlim(ax, datenum([xmin xmax]));
ax.XLimMode = 'manual';           % ← 新增：锁住范围，后续不被改写
xt = datenum(datetime(year(xmin),1,1):calmonths(2):xmax);
set(ax,'XTick',xt);
datetick(ax,'x','yyyy-mmm','keepticks','keeplimits');
ax.XTickMode = 'manual';          % ← 新增：锁住刻度
grid(ax,'on'); ax.GridAlpha=0.10; ax.MinorGridAlpha=0.08;
ax.YTick=[]; ax.Layer='top'; ax.FontSize=FS;
end



function draw_left_labels(ax, provOrder, blocks, blockNames, ...
                          rowStep, FS, leftMarginDays, labelShiftDays, xmin, ...
                          regIdxPerProv, regColors)
Nprov   = numel(provOrder);
x0      = datenum(xmin);
x_prov  = x0 - leftMarginDays;                  % 省名位置（靠近图区，右对齐）
x_reg   = x_prov - labelShiftDays;              % 区域名更靠左

for ip = 1:Nprov
    y0 = (Nprov - ip) * rowStep;
    c  = [0.12 0.12 0.12];
    if exist('regIdxPerProv','var') && ~isempty(regIdxPerProv)
        ii = regIdxPerProv(ip);
        if ii>=1 && ii<=size(regColors,1), c = regColors(ii,:); end
    end
    text(ax, x_prov, y0, char(provOrder(ip)), ...
        'HorizontalAlignment','right', 'VerticalAlignment','middle', ...
        'FontSize',FS, 'Color',c, 'FontName','Helvetica');
end

% 区域块标签本身也用对应颜色
if ~isempty(blocks) && ~isempty(blockNames) && numel(blocks)==numel(blockNames)
    ycum = 0;
    for b = 1:numel(blocks)
        n = numel(blocks{b});
        ycent = (Nprov - (ycum + n/2)) * rowStep;
        cc = regColors(min(b,size(regColors,1)),:);
        text(ax, x_reg, ycent, blockNames{b}, ...
            'Rotation',90, 'HorizontalAlignment','center', ...
            'FontSize',FS, 'Color',cc, 'FontName','Helvetica');
        ycum = ycum + n;
    end
end
end






function [d2,y2]=mask_range(d,y,xmin,xmax)
m=(d>=xmin)&(d<=xmax)&isfinite(y); d2=d(m); y2=y(m);
end



function draw_one(ax,d,yN,y0,colLine,colFill,a,lw)
x=datenum(d); y=y0+yN;
px=[x(:); flipud(x(:))]; py=[y0+zeros(size(y(:))); flipud(y(:))];
patch('XData',px,'YData',py,'FaceColor',colFill,'FaceAlpha',a,'EdgeColor','none','Parent',ax);
plot(ax,x,y,'Color',colLine,'LineWidth',lw);
end



function lw=choose_lw(pk,q50,q75,q90)
if pk>=q90, lw=1.4; elseif pk>=q75, lw=1.1; elseif pk>=q50, lw=0.9; else, lw=0.7; end
end



function vline(ax,x,hfrac,col,style,lw)
xl = xlim(ax); yl = ylim(ax); yspan = yl(2)-yl(1);
plot(ax, datenum([x x]), yl(1)+[0 hfrac]*yspan, style, 'Color',col,'LineWidth',lw);
xlim(ax, xl);   % 仍保持
end



function q=safe_quantile(x,p)
v=x(isfinite(x)); if isempty(v), q=Inf; else, q=quantile(v,p); end
end



function C = ridge_row_cmap_local(N)
a=[0.98 0.85 0.30]*0.8; b=[0.30 0.45 0.70]; t=linspace(0,1,max(1,N))';
C = a + (b - a).*t;
end
