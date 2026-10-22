function make_pinn_gallery_and_detail(varargin)
% 读取各省各波数据 → 调用 integrate_seir_odeint_style → 画图
% ① “总览图”：30省缩略图（每省叠所有波：观测、ODE拟合、f曲线提示）
% ② （可选）30个省级详图：每波两面板（拟合 & f）

%% ------------ 参数 ------------
ip = inputParser;
addParameter(ip,'rootDir','.',@(s)ischar(s)||isstring(s));
addParameter(ip,'featuresMat','features_aug.mat',@(s)ischar(s)||isstring(s));
addParameter(ip,'outName','Fig_PINN_Gallery',@(s)ischar(s)||isstring(s));
% 视觉
addParameter(ip,'FontName','Helvetica',@(s)ischar(s)||isstring(s));
addParameter(ip,'FS_title',10,@(x)isnumeric(x)&&isscalar(x));
addParameter(ip,'FS_axis',8,@(x)isnumeric(x)&&isscalar(x));
addParameter(ip,'FS_small',7.5,@(x)isnumeric(x)&&isscalar(x));
addParameter(ip,'LW_data',20,@(x)isnumeric(x)&&isscalar(x));
addParameter(ip,'LW_fit',1.2,@(x)isnumeric(x)&&isscalar(x));
addParameter(ip,'LW_f',0.9,@(x)isnumeric(x)&&isscalar(x));
addParameter(ip,'ShowFcurve',true,@islogical);
addParameter(ip,'MonthTick',1,@(x)isnumeric(x)&&isscalar(x)&&x>=1); % 1=每月, 2=隔月
addParameter(ip,'galleryWaveNo',1,@(x)isnumeric(x)&&isscalar(x)&&x>=1); % 总览图仅画这一波

addParameter(ip,'LegendOnOverview',true,@islogical);                  % 总览页是否放图例
addParameter(ip,'LegendModeProvince','page',@(s) any(strcmpi(s,{'page','axes','none'})));
addParameter(ip,'LegendLocation','southoutside',@(s)ischar(s)||isstring(s)); % 图例位置
addParameter(ip,'LegendFontSize',[],@(x) isempty(x) || (isnumeric(x)&&isscalar(x))); % 留空则用 FS_axis
addParameter(ip,'LegendToken',16,@(x)isnumeric(x)&&isscalar(x)&&x>0);          % 色块大小


addParameter(ip,'ShowFcurveRightAxis',true,@islogical);  % f 曲线放右轴 [0,1]
addParameter(ip,'ClipNegativeOde',true,@islogical);      % 绘图时将 I_ode<0 裁成 0
addParameter(ip,'OverviewContentType','image',@(s) any(strcmpi(s,{'image','vector'})));



% 功能
addParameter(ip,'smoothDays',1,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(ip,'exportPerProvince',false,@islogical);
parse(ip,varargin{:});
P = ip.Results;

rootDir = char(P.rootDir);
S = load(char(P.featuresMat));
T = get_features_table(S);
T.province = string(T.province);
[~, T.wave_no] = parse_wave_id(string(T.wave_id));

% 省序与颜色
provList = sort(unique(T.province,'stable'));
Nprov = numel(provList);
C = palette_nature30(Nprov);
prov2col = containers.Map(cellstr(provList), num2cell(C,2));

%% ------------ 预读所有波，组织到省缓存 ------------
Prov = struct();  % Prov(p).waves(k) 保存每个省每一波
for p = 1:Nprov
    prov = provList(p);
    Prov(p).name = prov;

    % 用“带字段的空结构体”初始化，避免第一次 end+1 时报错
    Prov(p).waves = struct( ...
        'wave_no',   {}, ...
        'date_obs',  {}, 'y_obs',  {}, ...
        'date_ode',  {}, 'y_ode',  {}, ...
        'f',         {}, 'date_f', {} );

    % 该省在 Features 中的所有行
    rows = find(T.province==prov);
    for r = rows'
        wn = T.wave_no(r);
        waveDir = fullfile(rootDir, char(prov), sprintf('Results%02d', wn));

        % ---- 读波数据：daily_cases、time_series(初值+f)、parameters(beta,gamma)
        [D, f, init, beta, gamma] = read_wave_inputs(waveDir, P.smoothDays);
        if isempty(D.date) || isempty(f) || ~isfinite(beta) || ~isfinite(gamma)
            warning('%s/Results%02d 缺少数据或参数（已跳过）。', prov, wn);
            continue;
        end

        % ---- 调用你的 ODE 求解器（先积分，再取 I_ode）
        out = integrate_seir_odeint_style(f, init, beta, gamma);

        % ---- 生成作图所需时间轴与序列
        date0 = D.date(1);
        d_ode = date0 + days(out.t(:));      % t=0..T-1 → 日历日期
        yode  = out.I_ode(:);
        if P.ClipNegativeOde
            yode = max(0, yode);             % 仅用于作图，避免负值撑爆坐标
        end
        df = date0 + days((0:numel(f)-1))';  % f 的日期轴

        % ---- 入缓存
        k = numel(Prov(p).waves) + 1;
        Prov(p).waves(k) = struct( ...
            'wave_no',  wn, ...
            'date_obs', D.date(:), 'y_obs',  D.y(:), ...
            'date_ode', d_ode(:),  'y_ode',  yode, ...
            'f',        f(:),      'date_f', df );
    end

    % ---- 按波号排序，确保后续使用的稳定性
    if ~isempty(Prov(p).waves)
        [~,ord] = sort([Prov(p).waves.wave_no]);
        Prov(p).waves = Prov(p).waves(ord);
    end
end


%% ------------ ① 30省总览图（每省叠所有波） ------------
fig = figure('Color','w','Name','PINN Gallery','Units','pixels','Position',[60 60 1600 1000]);
tlo = tiledlayout(fig, 5, 6, 'TileSpacing','compact','Padding','compact');

for p = 1:Nprov
    nexttile(tlo); ax = gca; hold(ax,'on'); box(ax,'off');
    set(ax,'FontName',P.FontName,'FontSize',P.FS_axis);
    col = prov2col(char(provList(p)));

    % 找该省的指定波（默认 wave_no==1）
    idx = find([Prov(p).waves.wave_no] == P.galleryWaveNo, 1);
    if isempty(idx)
        title(ax, sprintf('%s — Wave %02d (missing)', provList(p), P.galleryWaveNo), ...
            'FontSize',P.FS_axis+0.5,'FontWeight','normal');
        continue;
    end
    W = Prov(p).waves(idx);

    % ---- 左轴：观测 vs ODE（保证看得清）----
    yyaxis(ax,'left');
    plot(ax, W.date_obs, W.y_obs, '.', 'Color',[0.55 0.55 0.55], 'MarkerSize', P.LW_data);
    plot(ax, W.date_ode, W.y_ode, '-', 'Color', col, 'LineWidth', P.LW_fit);
    % 用分位数做稳健 y 轴，避免偶发尖峰“撑爆”
    yL = [W.y_obs; W.y_ode]; yL = yL(isfinite(yL));
    if ~isempty(yL)
        lo = max(0, prctile(yL, 0.0));   % 也可用 1 分位
        hi = prctile(yL, 99.5);          % 稍微留点顶端
        if hi<=lo, hi=lo+1; end
        ylim(ax, [lo, hi]);
    end
    ylabel(ax,'Daily cases');

    % ---- 右轴：f(t) 原尺度 [0,1] ----
    if P.ShowFcurveRightAxis && ~isempty(W.f)
        yyaxis(ax,'right');
        plot(ax, W.date_f, W.f, '--', 'Color', col*0.75+0.25, 'LineWidth', P.LW_f);
        % ylim(ax,[0 1]); 
        ax.YColor = col*0.6+0.4; ylabel(ax,'f(t)');
        yyaxis(ax,'left'); % 切回左轴，避免后续样式污染
    end

    % 轴与标题
    xlim_auto(ax, [datetime(2020,1,1) datetime(2022,12,31)]);
    % 设定刻度位置（每月或隔月）
    rngx = xlim(ax);
    if isa(rngx,'datetime') || isdatetime(rngx)
        xt = dateshift(rngx(1),'start','month'):calmonths(P.MonthTick):dateshift(rngx(2),'start','month');
        xticks(ax, xt);
    end
    % 直接用数字格式（无中文/英文名）
    xtickformat(ax,'yyyy-MM');      % 例如 2022-02
    ax.XTickLabelRotation = 30;

    rngx = xlim(ax);
    if isa(rngx,'datetime') || isdatetime(rngx)
        xt = dateshift(rngx(1),'start','month'):calmonths(P.MonthTick):dateshift(rngx(2),'start','month');
        xticks(ax, xt);
    end
    ax.XTickLabelRotation=30; grid(ax,'on'); ax.GridAlpha=0.10;

    title(ax, sprintf('%s — Wave %02d', provList(p), P.galleryWaveNo), ...
        'FontSize',P.FS_axis+0.5,'FontWeight','normal','Color',[0.12 0.12 0.12]);


end
title(tlo, "PINN fits (observed vs ODE) with scaled f-curves — by province", ...
    'FontSize', P.FS_title, 'FontWeight','bold');

% 导出总览
savefig(fig, [char(P.outName) '.fig']);
exportgraphics(fig, [char(P.outName) '.pdf'], 'ContentType', P.OverviewContentType,'BackgroundColor','white');

%% ========== A4 画册式 PDF：第1页总览 + 后续30页每省 ==========

outPdf = [char(P.outName) '_Booklet.pdf'];
export_booklet_A4(outPdf, Prov, provList, prov2col, P);


function export_booklet_A4(outPdf, Prov, provList, prov2col, P)

% ---------- 纸张和字体 ----------
A4px = [0 0 2480 3508];   % 近似像素尺寸（300dpi）；只是做个宽高比
% 也可以用 Position 设一个接近 A4 比例的像素窗口，exportgraphics 会自适应

% ---------- 第 1 页：30省总览（6×5，第一波） ----------
% ---------- 第 1+ 页：总览（每页 6×5，只画指定的第一波） ----------
rows = 6; cols = 5; perPage = rows*cols;
Nprov = numel(provList);
NgalPages = ceil(Nprov / perPage);

for g = 1:NgalPages
    idxStart = (g-1)*perPage + 1;
    idxEnd   = min(g*perPage, Nprov);

    figG = figure('Color','w','Units','pixels','Position',[100 50 1240 1754], ...
                  'Name',sprintf('Gallery p.%d', g), 'Visible','off');
    tloG = tiledlayout(figG, rows, cols, 'TileSpacing','compact','Padding','compact');
    axList = gobjects(0);

    for p = idxStart:idxEnd
        nexttile(tloG); ax = gca; axList(end+1) = ax; %#ok<AGROW>
        hold(ax,'on'); box(ax,'off');
        set(ax,'FontName',P.FontName,'FontSize',P.FS_axis);
        col = prov2col(char(provList(p)));

        % 该省的指定波（默认 Wave 01）
        iW = find([Prov(p).waves.wave_no]==P.galleryWaveNo,1);
        if isempty(iW)
            title(ax, sprintf('%s — Wave %02d (NA)', provList(p), P.galleryWaveNo), ...
                'FontSize',P.FS_axis); axis(ax,'off'); continue;
        end
        W = Prov(p).waves(iW);

        % 左轴：观测 vs ODE
        yyaxis(ax,'left');
        plot(ax, W.date_obs, W.y_obs, '.', 'Color',[0.52 0.52 0.52], 'Markersize', P.LW_data);
        plot(ax, W.date_ode, W.y_ode, '-', 'Color', col, 'LineWidth', P.LW_fit);
        yL = [W.y_obs; W.y_ode]; yL = yL(isfinite(yL));
        if ~isempty(yL)
            lo = max(0, prctile(yL,0.0)); hi = prctile(yL,99.5);
            if hi<=lo, hi=lo+1; end, ylim(ax,[lo,hi]);
        end
        ylabel(ax,'Daily cases');

        % 右轴：f(t)
        if P.ShowFcurveRightAxis && ~isempty(W.f)
            yyaxis(ax,'right');
            plot(ax, W.date_f, W.f, '--', 'Color', col*0.75+0.25, 'LineWidth', P.LW_f);
            % ylim(ax,[0 1]); 
            ax.YColor = col*0.6+0.4; ylabel(ax,'f(t)');
            yyaxis(ax,'left');
        end

        % 时间轴（用纯数字避免“#”字形问题）
        xlim_auto(ax,[datetime(2020,1,1) datetime(2022,12,31)]);
        xtickformat(ax,'yyyy-MM');
        rngx = xlim(ax);
        if isa(rngx,'datetime') || isdatetime(rngx)
            xt = dateshift(rngx(1),'start','month'):calmonths(P.MonthTick):dateshift(rngx(2),'start','month');
            xticks(ax, xt);
        end
        ax.XTickLabelRotation=30; grid(ax,'on'); ax.GridAlpha=0.10;

        % 标题加“跳转页码提示”（总览第 g 页 → 省级页序号：p + 1）
        pageTarget = p + 1;  % 省级页从第2页开始
        title(ax, sprintf('%s — Wave %02d  (p.%02d)', provList(p), P.galleryWaveNo, pageTarget), ...
            'FontSize',P.FS_axis+0.3,'FontWeight','normal','Color',[0.15 0.15 0.15]);
    end

    title(tloG, "PINN fits (observed vs ODE, wave 01) with f(t) on right axis — by province", ...
          'FontSize', P.FS_title, 'FontWeight','bold');

    % —— 本页图例（只放一份）：挂在第一个轴，再放到页脚，不占格子
    if P.LegendOnOverview && ~isempty(axList)
        ax0 = axList(1);
        hold(ax0,'on');
        hData = plot(ax0, NaN,NaN, '.', 'Color',[0.52 0.52 0.52], 'Markersize', P.LW_data, 'DisplayName','Real data');
        hFit  = plot(ax0, NaN,NaN, '-', 'Color',[0.20 0.20 0.20], 'LineWidth', P.LW_fit,  'DisplayName','Fitting');
        hF    = plot(ax0, NaN,NaN, '--', 'Color',[0.35 0.35 0.35], 'LineWidth', P.LW_f,    'DisplayName','f(t)');

        lgd_fs = ternary(isempty(P.LegendFontSize), P.FS_axis, P.LegendFontSize);
        lgdG = legend(ax0, [hData hFit hF], {'Real data','Fitting','f(t)'}, ...
            'Orientation','horizontal', 'Box','off', ...
            'FontName', P.FontName, 'FontSize', lgd_fs);
        lgdG.ItemTokenSize = [P.LegendToken P.LegendToken];
        lgdG.Layout.Tile   = 'south';   % 放到整页页脚
    end

    % —— 导出（总览页建议 image，加速；如果你坚持矢量，改成 'vector'）
    exportgraphics(figG, outPdf, 'ContentType', P.OverviewContentType,'BackgroundColor','white', ...
                   'Append', g>1);   % 第一页 Append=false，其余 true
    drawnow;
end

title(tloG, "PINN fits (observed vs ODE, wave 01) with f(t) on right axis", ...
      'FontSize', P.FS_title, 'FontWeight','bold');

% 导出第1页（矢量）
exportgraphics(figG, outPdf, 'ContentType','vector','BackgroundColor','white');
drawnow;

% ---------- 第 2～31 页：每省 1 页，自适应网格 ----------
for p = 1:numel(provList)
    waves = Prov(p).waves;
    nW = numel(waves);
    if nW==0
        figP = figure('Color','w','Units','pixels','Position',[120 60 1240 1754], ...
                      'Name',sprintf('%s — NA',provList(p)), 'Visible','off');
        tloP = tiledlayout(figP,1,1,'TileSpacing','compact','Padding','compact');
        nexttile(tloP); axis off;
        text(0.5,0.5,sprintf('%s (no waves)',provList(p)), ...
             'Units','normalized','HorizontalAlignment','center','FontSize',P.FS_title);
        exportgraphics(figP, outPdf, 'ContentType','vector','BackgroundColor','white','Append',true);
        continue;
    end

    [R,C] = grid_for_nwaves(nW);
    figP = figure('Color','w','Units','pixels','Position',[120 60 1240 1754], ...
                  'Name',sprintf('%s — %d waves',provList(p),nW), 'Visible','off');
    tloP = tiledlayout(figP, R, C, 'TileSpacing','compact','Padding','compact');

    col = prov2col(char(provList(p)));
    axList = gobjects(0);   % <<< 初始化收集本页的轴

    for k = 1:nW
        nexttile(tloP);
        ax = gca; axList(end+1) = ax;   % <<< 收集
        hold(ax,'on'); box(ax,'off');
        set(ax,'FontName',P.FontName,'FontSize',P.FS_axis);

        W = waves(k);

        % 左轴：观测 vs ODE
        yyaxis(ax,'left');
        plot(ax, W.date_obs, W.y_obs, '.', 'Color',[0.52 0.52 0.52], 'MarkerSize', P.LW_data);  % 代表点
        plot(ax, W.date_ode, W.y_ode, '-', 'Color', col, 'LineWidth', P.LW_fit);
        yL = [W.y_obs; W.y_ode]; yL = yL(isfinite(yL));
        if ~isempty(yL)
            lo = max(0, prctile(yL, 0.0)); hi = prctile(yL, 99.5);
            if hi<=lo, hi=lo+1; end, ylim(ax,[lo,hi]);
        end
        ylabel(ax,'Daily cases');

        % 右轴：f(t)
        if P.ShowFcurveRightAxis && ~isempty(W.f)
            yyaxis(ax,'right');
            plot(ax, W.date_f, W.f, '--', 'Color', col*0.75+0.25, 'LineWidth', P.LW_f);
            ax.YColor = col*0.6+0.4; ylabel(ax,'f(t)');
            yyaxis(ax,'left');
        end

        xlim_auto(ax,[datetime(2020,1,1) datetime(2022,12,31)]);
        rngx = xlim(ax);
        if isa(rngx,'datetime') || isdatetime(rngx)
            xt = dateshift(rngx(1),'start','month'):calmonths(P.MonthTick):dateshift(rngx(2),'start','month');
            xticks(ax, xt);
        end
        xtickformat(ax,'yyyy-MM');
        ax.XTickLabelRotation = 30;

        grid(ax,'on'); ax.GridAlpha=0.10;
        title(ax, sprintf('%s — Wave %02d', provList(p), W.wave_no), ...
              'FontSize',P.FS_axis+0.3,'FontWeight','normal');
    end

    % —— 省级页图例（每页一份，放页脚；注意：放在循环之后！） ——
    if ~strcmpi(P.LegendModeProvince,'none') && ~isempty(axList)
        ax0 = axList(1);  hold(ax0,'on');  % 宿主轴，不占 tile
        hData = plot(ax0, NaN,NaN, '.', 'Color',[0.52 0.52 0.52], 'MarkerSize', P.LW_data, 'DisplayName','Real data');
        hFit  = plot(ax0, NaN,NaN, '-', 'Color', col,                 'LineWidth', P.LW_fit, 'DisplayName','Fitting');
        hF    = plot(ax0, NaN,NaN, '--','Color', col*0.75+0.25,       'LineWidth', P.LW_f,   'DisplayName','f(t)');

        lgd_fs = ternary(isempty(P.LegendFontSize), P.FS_axis, P.LegendFontSize);
        lgdP = legend(ax0, [hData hFit hF], {'Real data','Fitting','f(t)'}, ...
            'Orientation','horizontal', 'Box','off', ...
            'FontName', P.FontName, 'FontSize', lgd_fs);
        lgdP.ItemTokenSize = [P.LegendToken P.LegendToken];
        lgdP.Layout.Tile   = 'south';   % 放到整页页脚
    end

    % 页标题与导出
    title(tloP, sprintf('%s — PINN fits (observed vs ODE) with f(t)  |  n=%d waves', ...
          provList(p), nW), 'FontSize', P.FS_title, 'FontWeight','bold');

    exportgraphics(figP, outPdf, 'ContentType','vector','BackgroundColor','white','Append',true);
end


%     % 省级页标题
%     title(tloP, sprintf('%s — PINN fits (observed vs ODE) with f(t) (n=%d waves)', provList(p), nW), ...
%         'FontSize', P.FS_title, 'FontWeight','bold');
% 
%     % 追加导出到同一 PDF
%     exportgraphics(figP, outPdf, 'ContentType','vector','BackgroundColor','white','Append',true);
%     close(figP);
% end
end

% ---- 行列规则（≤4 →2x2；5–6→3x2；7–9→3x3；10–16→5x3）----
function [R,C] = grid_for_nwaves(n)
if n<=4, R=2; C=2;
elseif n<=6, R=3; C=2;
elseif n<=9, R=3; C=3;
else, R=5; C=3;  % up to 15
end
end



end

%% =================== I/O 与工具函数 ===================
function T = get_features_table(S)
T = [];
if isfield(S,'Features'), T = S.Features; else
    f = fieldnames(S);
    for i=1:numel(f), if istable(S.(f{i})), T = S.(f{i}); break; end, end
end
assert(~isempty(T),'features_aug.mat 中未找到表。');
% 关键列名容错
vn = lower(string(T.Properties.VariableNames));
req = {'province','wave_id'};
for k=1:numel(req)
    assert(any(vn==req{k}), 'Features 缺少列：%s', req{k});
end
if ~any(vn=="start_date"), T.start_date = NaT(height(T),1); end %#ok<AGROW>
if ~isdatetime(T.start_date), try, T.start_date = datetime(T.start_date); catch, end, end
end

function [D, f, init, beta, gamma] = read_wave_inputs(waveDir, smoothDays)
% 更鲁棒地读取：
% - daily_cases.xlsx（日期+日增）
% - results/ 或 Results/ 下的 time_series*.xlsx（初值 + f(t)）
% - results/ 或 Results/ 下的 parameters*.xlsx（beta, gamma）
%
% 列名可能为中文/带符号/希腊字母，本函数会做“表头归一化”并多策略匹配。
%
D = struct('date',datetime.empty,'y',[]);
f = []; init = struct('Se0',NaN,'S0',NaN,'I0',NaN,'R0',NaN);
beta = NaN; gamma = NaN;

%% ---- 1) daily_cases.xlsx ----
dc = fullfile(waveDir,'daily_cases.xlsx');
if isfile(dc)
    T = readtable(dc,'VariableNamingRule','preserve');
    if size(T,1) >= 2
        d = T{2:end,1};
        if ~isdatetime(d)
            % 多尝试几种常见格式
            try, d = datetime(d);
            catch, try, d = datetime(d,'InputFormat','yyyy-M-d');
            catch, d = datetime(d,'InputFormat','yyyy-MM-dd'); end, end
        end
        y = double(T{2:end,2});
        if smoothDays>1
            y = movmean(y,[floor((smoothDays-1)/2) ceil((smoothDays-1)/2)],'omitnan');
        end
        D.date = d(:); D.y = y(:);
    end
end

%% ---- 2) time_series*.xlsx：results/Results 目录下模糊查找 ----
candDirs = {fullfile(waveDir,'results'), fullfile(waveDir,'Results'), waveDir};
tsFile = '';
for dd = 1:numel(candDirs)
    if ~isfolder(candDirs{dd}), continue; end
    L = dir(fullfile(candDirs{dd}, '*time*series*.xlsx'));
    if isempty(L), L = dir(fullfile(candDirs{dd}, 'time_series.xlsx')); end
    if ~isempty(L), [~,i] = max([L.datenum]); tsFile = fullfile(L(i).folder, L(i).name); break; end
end

if ~isempty(tsFile) && isfile(tsFile)
    TT = readtable(tsFile,'VariableNamingRule','preserve');
    % 归一化表头 → 更易匹配
    rawNames = string(TT.Properties.VariableNames);
    norm = normalize_names(rawNames);

    % 尝试匹配 Se/S/I/R/f 的列
    Se_col = find_first(norm, ["se","se_0","se0","e","expo","exposed"]);
    S_col  = find_first(norm, ["s","s_0","s0","sus","susceptible"]);
    I_col  = find_first(norm, ["i","i_0","i0","inf","infectious"]);
    R_col  = find_first(norm, ["r","r_0","r0","rec","removed","recovered"]);
    f_col  = find_first(norm, ["f","ft","f_t","phi","policy","u","control"]);

    % 初值取第 1 行
    if ~isempty(Se_col), init.Se0 = safe_num(TT{1,Se_col}); end
    if ~isempty(S_col),  init.S0  = safe_num(TT{1,S_col});  end
    if ~isempty(I_col),  init.I0  = safe_num(TT{1,I_col});  end
    if ~isempty(R_col),  init.R0  = safe_num(TT{1,R_col});  end

    % f 序列：若没匹配到，启用“候选列扫描（0~1 区间的长数值列）”
    if ~isempty(f_col)
        f = double(TT{:,f_col});
    else
        % 从所有数值列中挑选“在 [0,1.2] 区间、长度>=5”的列
        numCols = find(varfun(@isnumeric,TT,'OutputFormat','uniform'));
        best = [];
        for c = numCols
            vv = double(TT{:,c});
            if nnz(isfinite(vv)) >= 5
                lo = prctile(vv,1,'all','omitnan'); hi = prctile(vv,99,'all','omitnan');
                if lo >= -0.05 && hi <= 1.20
                    best = c; break;
                end
            end
        end
        if ~isempty(best), f = double(TT{:,best}); end
    end
end

%% ---- 3) parameters*.xlsx（beta/gamma）----
pmFile = '';
for dd = 1:numel(candDirs)
    if ~isfolder(candDirs{dd}), continue; end
    L = dir(fullfile(candDirs{dd}, '*param*.xlsx'));
    if isempty(L), L = dir(fullfile(candDirs{dd}, 'parameters.xlsx')); end
    if ~isempty(L), [~,i] = max([L.datenum]); pmFile = fullfile(L(i).folder, L(i).name); break; end
end

if ~isempty(pmFile) && isfile(pmFile)
    PP = readtable(pmFile,'VariableNamingRule','preserve');

    % 两种结构都支持：
    %  A) 宽表：有名为 beta/gamma 的列（或变体）
    %  B) 列表：第一列是 Name，第二列是 Value

    % ---- 情形 A：宽表 ----
    rawP = string(PP.Properties.VariableNames);
    normP = normalize_names(rawP);

    bcol = find_first(normP, ["beta","b","β","beta_","par_beta"]);
    gcol = find_first(normP, ["gamma","g","γ","gamma_","par_gamma"]);

    if ~isempty(bcol) && height(PP)>=1, beta  = safe_num(PP{1,bcol}); end
    if ~isempty(gcol) && height(PP)>=1, gamma = safe_num(PP{1,gcol}); end

    % ---- 情形 B：Name/Value 列表（如果上面没读到）----
    if (~isfinite(beta) || ~isfinite(gamma)) && width(PP)>=2
        names = normalize_text_col(PP{:,1});
        vals  = PP{:,2};
        if ~isfinite(beta)
            idx = find(contains(names, ["beta","b","β"]),1);
            if ~isempty(idx), beta = safe_num(vals(idx)); end
        end
        if ~isfinite(gamma)
            idx = find(contains(names, ["gamma","g","γ"]),1);
            if ~isempty(idx), gamma = safe_num(vals(idx)); end
        end
    end
end
end

% ====== 小工具 ======
function n = safe_num(x)
if iscell(x), x = x{1}; end
if isstring(x) || ischar(x), x = str2double(x); end
n = double(x);
end

function nn = normalize_names(vnames)
% 小写 → 去空格/括号/单位/点号/百分号等 → latin 化希腊字母
nn = lower(string(vnames));
nn = regexprep(nn,'[\s%()\[\]\{\}\-\.\\/]+','');   % 去分隔符
% 希腊字母替换
nn = strrep(nn,'β','beta'); nn = strrep(nn,'γ','gamma'); nn = strrep(nn,'φ','phi');
% 常见别名统一
nn = strrep(nn,'susceptibles','susceptible');
nn = strrep(nn,'infected','infectious');
nn = strrep(nn,'removed','recovered');
end

function s = normalize_text_col(col)
% 把第一列名称归一化（用于 Name/Value 列表）
if iscell(col) || isstring(col) || ischar(col)
    s = lower(string(col));
else
    s = lower(string(string(col)));
end
s = regexprep(s,'[\s%()\[\]\{\}\-\.\\/]+','');
s = strrep(s,'β','beta'); s = strrep(s,'γ','gamma'); s = strrep(s,'φ','phi');
end

function idx = find_first(names, keys)
idx = [];
for k=1:numel(keys)
    j = find(names == keys(k), 1);
    if ~isempty(j), idx=j; return; end
end
for k=1:numel(keys)
    j = find(contains(names, keys(k)), 1);
    if ~isempty(j), idx=j; return; end
end
end


function idx = pick_col(cols, keys)
idx = [];
for k=1:numel(keys)
    j = find(cols==keys(k),1);
    if ~isempty(j), idx=j; return; end
end
% 次优：contains
for k=1:numel(keys)
    j = find(contains(cols, keys(k)),1);
    if ~isempty(j), idx=j; return; end
end
end

function xlim_auto(ax, defaultRange)
if nargin<2 || isempty(defaultRange)
    defaultRange = [datetime(2020,1,1) datetime(2022,12,31)];
end
xl = xlim(ax);

% datetime 轴
if isa(xl,'datetime') || isdatetime(xl)
    if any(isnat(xl)) || diff(xl)<=duration(0,0,0)
        xlim(ax, defaultRange);
    end
% 数值轴（datenum）
else
    if any(~isfinite(xl)) || xl(1)>=xl(2)
        xlim(ax, datenum(defaultRange));
    end
end
end


function [prov,no]=parse_wave_id(wid)
wid=string(wid); prov=strings(numel(wid),1); no=zeros(numel(wid),1);
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

function C = palette_nature30(N)
% 自然风 30 色
if nargin<1, N=30; end
C=zeros(N,3); phi=137.508; h0=15;
for i=1:N
    h = mod(h0+(i-1)*phi,360);
    s = 0.55 + 0.10*sin((i-1)*pi/6);
    v = 0.78 + 0.08*cos((i-1)*pi/8);
    C(i,:) = hsv2rgb([h/360 s v]);
end
C = min(max(C,0),1);
end


function y = ternary(cond, a, b)
if cond, y = a; else, y = b; end
end
