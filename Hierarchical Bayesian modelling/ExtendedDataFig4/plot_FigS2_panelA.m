function H = plot_FigS2_panelA(varargin)
% Panel S2(a) — Province radar (inner=province RE; outer=variant-adjusted by weights)
% 用法示例：
% H = plot_FigS2_panelA('results_root','results_hbm', ...
%        'features_mat','features_aug.mat', ...
%        'target','delta_early', ...
%        'weightScheme','waves', ...       % waves|cases|duration
%        'radarMaxScale',[], ...           % 自动nice刻度；可手动如0.6
%        'radialTransform','asinh', ...    % 'identity'|'sqrt'|'asinh'
%        'softness',0.25, ...              % asinh的软阈
%        'labelEvery',18, ...              % 约每 N/labelEvery 标1个省名；设1为全部
%        'renderer','painters');           % 'painters'(推荐) | 'opengl'

%% ---------- 参数 ----------
p = inputParser;
p.addParameter('results_root','results_hbm',@(s)ischar(s)||isstring(s));
p.addParameter('features_mat','features_aug.mat',@(s)ischar(s)||isstring(s));
p.addParameter('target','',@(s)ischar(s)||isstring(s));
p.addParameter('weightScheme','waves',@(s)ischar(s)||isstring(s)); % waves|cases|duration
p.addParameter('radarMaxScale',[],@(x)isnumeric(x)||isempty(x));
p.addParameter('radialTransform','identity',@(s)ischar(s)||isstring(s)); % 'identity'|'sqrt'|'asinh'
p.addParameter('softness',0.25,@(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('labelEvery',18,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('renderer','painters',@(s)ischar(s)||isstring(s));

p.addParameter('labelRFactor', 1.0, @(x)isnumeric(x)&&isscalar(x)&&x>0);   % 省名半径相对 r_plot_max 的倍数，默认 1.06
p.addParameter('labelTangentOffset', 0.002, @(x)isnumeric(x)&&isscalar(x)&&x>=0); % 切向微偏移（避免贴边）

p.parse(varargin{:});
args = p.Results;

%% ---------- 颜色 ----------
Cgrp  = safe_palette(@() smf_colorlib('groups8'),   lines(8));         % region配色
CvarO = [0.05 0.35 0.80];                                              % 外环均值点（蓝）
Cgray = [0.20 0.20 0.20];
fontMain = pick_font({'Arial','Helvetica','Nimbus Sans L'});           % 清晰、不发粗

%% ---------- 载入 cache / Features ----------
cache = load_cache(args.results_root);
prov_levels = string(cache.prov_names(:));
targets_all = string(cache.active_targets(:)');

target = select_target(args.target, targets_all);

T = load_features(args.features_mat);
[colProv,colVariant,colRegion,colCases,colDuration] = resolve_cols(T);

% 标准化 variant、region
T.variant = canon_variant(string(T.(colVariant)));
T.variant = replace(T.variant,"Original","Alpha"); % 与 ranef 对齐
region_series = region_series_of(T, colRegion);
region_of_prov = province_region_mode(T, prov_levels, region_series);

%% ---------- 读取 ranef ----------
[mu_prov, sd_prov, ciP, mu_var, sd_var, var_levels] = read_ranef(args.results_root, target, prov_levels);

%% ---------- 权重 W_{i,k} ----------
W = compute_weights(T, prov_levels, var_levels, lower(args.weightScheme), colCases, colDuration); % P×K

%% ---------- 严格尺寸校验 ----------
mu_prov = double(mu_prov(:));         % P×1
sd_prov = double(sd_prov(:));         % P×1
mu_var  = double(mu_var(:));          % K×1
sd_var  = double(sd_var(:));          % K×1
[Pw, Kw] = size(W);
assert(Pw==numel(mu_prov), 'W 的行数(P=%d) 与 mu_prov 长度(%d) 不一致。', Pw, numel(mu_prov));
assert(Kw==numel(mu_var),  'W 的列数(K=%d) 与 mu_var 长度(%d) 不一致。', Kw, numel(mu_var));

%% ---------- 外环：调整后均值与方差传播 ----------
adj_variant = W * mu_var;                         % P×1
adj_sd      = sqrt( (W.^2) * (sd_var.^2) );       % P×1
mu_adj = mu_prov + adj_variant;                   % P×1
sd_adj = sqrt( sd_prov.^2 + adj_sd.^2 );          % P×1

%% ---------- 排序（按 region 分块，块内按 mu_prov 降序） ----------
[regBlocks, order_idx, reg_names] = block_order_by_region(region_of_prov, mu_prov);

%% ---------- rmax（原始尺度） ----------
rmax_auto = max(abs([mu_prov; mu_adj; ciP(:)]), [], 'omitnan');
rmax = args.radarMaxScale; 
if isempty(rmax) || ~isfinite(rmax), rmax = nice_rmax(rmax_auto); end

%% ---------- 非线性径向映射（缓解中心拥挤；不改符号/顺序） ----------
tf = lower(string(args.radialTransform));
s  = args.softness;          % asinh 的“软阈”
RMAX = rmax;                 % 原始最大半径（用于归一化）
switch tf
    case "sqrt"
        rmap = @(x) sign(x).*sqrt(abs(x)./max(RMAX,eps)) .* RMAX;
        labelTF = 'sqrt radial scale';
    case "asinh"
        denom = asinh(RMAX/max(s,1e-9));
        rmap = @(x) sign(x).* (asinh(abs(x)./max(s,1e-9))./max(denom,1e-9)) .* RMAX;
        labelTF = sprintf('asinh radial scale (s=%.2g)', s);
    otherwise
        rmap = @(x) x;
        labelTF = 'linear radial scale';
end

%% ---------- 绘图 ----------
fig = figure('Color','w','Name','FigS2(a) Province radar');
set(fig,'Renderer',char(args.renderer));   % 'painters' 推荐出版
ax  = axes(fig); hold(ax,'on'); axis(ax,'off');
ax.SortMethod = 'childorder';              % 保持绘制顺序（对 opengl 更稳定）

% title(ax, sprintf('Province radar'), ...
%     'FontName',fontMain,'FontSize',15,'FontWeight','normal','Color',Cgray);
% —— 可选中的标题（annotation 版）——
% ttl = annotation(fig, 'textbox', [0.10 0.93 0.80 0.05], ...  % [x y w h] in figure normalized units
%     'String','Province radar', ...
%     'EdgeColor','none', 'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
%     'FontName', fontMain, 'FontSize', 15, 'FontWeight','normal', 'Color', Cgray, ...
%     'Interpreter','none', 'HitTest','on', 'PickableParts','all');


N = numel(order_idx);
theta  = linspace(0,2*pi,max(N,1)+1)'; theta(end) = [];
dtheta = 2*pi/max(1,N);

% 网格圈（按映射后的半径绘制）
for rr=1:4
    draw_circle(ax, rmap(rr/4*rmax), [0.85 0.85 0.85]);
end

% 扇区
for kk=1:N
    i = order_idx(kk);
    ang1 = theta(kk)-dtheta/2; ang2 = theta(kk)+dtheta/2;
    colReg = region_color(Cgrp, regBlocks(kk), numel(reg_names));

    % 内环 CI 楔形（映射）
    if all(isfinite(ciP(i,:)))
        draw_wedge(ax, ang1, ang2, rmap(ciP(i,1)), rmap(ciP(i,2)), colReg, 0.22);
    end
    % 内环均值点
    draw_mark(ax, (ang1+ang2)/2, rmap(mu_prov(i)), 'k');

    % 外环：均值 + 近似CI（映射）
    r1 = mu_adj(i); r1s = sd_adj(i);
    if isfinite(r1s) && r1s>0
        draw_wedge(ax, ang1, ang2, rmap(r1-1.96*r1s), rmap(r1+1.96*r1s), [0.2 0.2 0.2], 0.10);
    end
    draw_mark(ax, (ang1+ang2)/2, rmap(r1), CvarO);
end

% 边界与分块分隔
r_plot_max = max(abs([rmap(-rmax) rmap(rmax)]));
pad = 1.4;
axis(ax,'equal'); xlim(ax,[-pad pad]*r_plot_max); ylim(ax,[-pad pad]*r_plot_max);
draw_separators(ax, theta, dtheta, regBlocks, r_plot_max);

%% —— 省名（左右半区、可全部显示、强制非加粗 + 清晰字体）——
step  = max(1, round(N / max(1,args.labelEvery)));  % labelEvery=1 => step=1 => 全部显示
step = 1;
for kk = 1:N
    % 每隔 step 个标一次；step=1 => 全部
    if mod(kk-1, step) ~= 0, continue; end

    i = order_idx(kk);
    ang1 = theta(kk)-dtheta/2; 
    ang2 = theta(kk)+dtheta/2;
    ang  = (ang1+ang2)/2;

    ang_n = atan2(sin(ang), cos(ang));  % (-pi,pi]
    rlab = args.labelRFactor * r_plot_max;
    [xlab, ylab] = pol2cart(ang_n, rlab);

    isRight = abs(ang_n) <= pi/2;
    if isRight
        rotDeg = rad2deg(ang_n);  hAlign = 'left';
    else
        rotDeg = rad2deg(ang_n) + 180;  hAlign = 'right';
    end

    off  = args.labelTangentOffset * rlab;  tnx = -sin(ang_n);  tny = cos(ang_n);
    xlab = xlab + (isRight*2-1)*off*tnx;
    ylab = ylab + (isRight*2-1)*off*tny;

    text(ax, xlab, ylab, char(prov_levels(i)), ...
        'FontName', fontMain, 'FontSize', 14, 'FontWeight','normal', 'FontAngle','normal', ...
        'Color',[0.2 0.2 0.2], ...
        'HorizontalAlignment',hAlign, 'VerticalAlignment','middle', ...
        'Rotation',rotDeg, 'Clipping','off', 'Interpreter','none');
end

%% 底部说明（annotation：可选中/拖动；不影响透明叠层）
noteStr = sprintf([ ...
    'Inner: province RE (mean, CI wedge)  |  Outer: variant-adjusted (mean, approx. CI)\n', ...
    '%s'], labelTF);   % 第二行：径向尺度说明

ann = annotation(fig,'textbox',[0.20 0.02 0.60 0.06], ...  % ← 缩窄为 60% 宽，居中
    'String', noteStr, 'EdgeColor','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontName',fontMain,'FontSize',11,'Color',[0.2 0.2 0.2], ...
    'Interpreter','none','FitBoxToText','off', ...    % 固定盒宽，启用换行
    'HitTest','on','PickableParts','all');


%% 图例（示例句柄）
hInnerDot = plot(ax, NaN, NaN, 'o', 'MarkerSize',4, 'MarkerFaceColor','k', ...
    'MarkerEdgeColor',[0.1 0.1 0.1], 'DisplayName','Inner: province mean');
hOuterDot = plot(ax, NaN, NaN, 'o', 'MarkerSize',4, 'MarkerFaceColor',CvarO, ...
    'MarkerEdgeColor',[0.1 0.1 0.1], 'DisplayName','Outer: adjusted mean');
hInnerBand = patch(ax, NaN, NaN, [0.3 0.3 0.3], 'FaceAlpha',0.42, 'EdgeColor','none', ...
    'DisplayName','Inner 95% CI (region color)');
hOuterBand = patch(ax, NaN, NaN, [0.2 0.2 0.2], 'FaceAlpha',0.30, 'EdgeColor','none', ...
    'DisplayName','Outer 95% CI (approx.)');
L = legend(ax, [hInnerDot hOuterDot hInnerBand hOuterBand], ...
    'Location','southoutside', 'Orientation','horizontal', ...
    'Box','off', 'FontName',fontMain,'FontSize',12, 'Interpreter','none');
L.NumColumns = 2;                 % ← 两列 => 两行（4项时每行2项）
L.ItemTokenSize = [12 12];        % ← 符号更紧凑
drawnow;                          % 让 legend 定位后再微调位置（可选）
% —— 如需更窄，可手动稍微减小宽度（可选）——
posL = L.Position;                % [x y w h] in normalized units
posL(3) = posL(3)*0.85;           % 宽度收窄 15%
posL(1) = 0.5 - posL(3)/2;        % 居中
L.Position = posL;


%% 返回句柄
H = struct('fig',fig,'ax',ax,'order_idx',order_idx,'region_names',reg_names);

%% ================= helpers（本文件内） =================

function cache = load_cache(root)
    S = load(fullfile(root,'cache','hbm_cache.mat')); cache=S.cache;
end

function T = load_features(fpath)
% 鲁棒 Features 读取器：
% - 先匹配常见名：Features/features/T/Table/tbl/Features_aug/features_aug
% - 若没匹配到，找“类型=table”的变量
% - 若仍没有，尝试结构体内嵌的 table 字段
% - 全失败时报错并列出变量名与类型
    if exist(fpath,'file')~=2
        error('未找到 features_mat 文件：%s', fpath);
    end
    info = whos('-file', fpath);
    if isempty(info), error('文件 %s 为空。', fpath); end
    cand = ["Features","features","T","Table","tbl","Features_aug","features_aug"];
    names = string({info.name}); lowerNames = lower(names);
    hit = [];
    for ii = 1:numel(cand)
        j = find(lowerNames == lower(cand(ii)), 1);
        if ~isempty(j), hit = names(j); break; end
    end
    if ~isempty(hit)
        S = load(fpath, char(hit));
        T0 = S.(char(hit));
        if istable(T0), T = T0; return; end
    end
    isTab = strcmp({info.class}, 'table');
    if any(isTab)
        pickName = info(find(isTab,1)).name;
        S = load(fpath, pickName);
        T = S.(pickName); return;
    end
    S = load(fpath);
    fn = fieldnames(S);
    for k = 1:numel(fn)
        val = S.(fn{k});
        if istable(val), T = val; return; end
        if isstruct(val)
            subfn = fieldnames(val);
            for m = 1:numel(subfn)
                vv = val.(subfn{m});
                if istable(vv), T = vv; return; end
            end
        end
    end
    listing = join( compose("%s : %s", names, string({info.class})), ", " );
    error('Features 缺失：%s 中未找到 table 变量。文件包含：%s', fpath, listing);
end

function [cp,cv,cr,cc,cd] = resolve_cols(T)
    cp = pick(T,["province","Province"]);
    cv = pick(T,["variant","Variant"]);
    cr = pick(T,["region category","region","Region"]);
    cc = pick(T,["cases","cases_total","cum_cases","n_cases","case_count"]);
    cd = pick(T,["duration","duration_days","len","length","Tlen","n_days"]);
    if strlength(cp)==0 || strlength(cv)==0
        error('Features 缺少 province / variant 列。');
    end
end

function s = pick(T,cand)
    vn = string(T.Properties.VariableNames);
    s = ""; 
    for ii=1:numel(cand)
        h = find(strcmpi(vn, cand(ii)), 1);
        if ~isempty(h), s = vn(h); return; end
    end
end

function target = select_target(t, pool)
    if strlength(t)==0, target=pool(1); else
        k0 = lower(regexprep(string(t),'[^a-z0-9]+',''));
        pp = lower(regexprep(pool,'[^a-z0-9]+',''));
        h = find(pp==k0,1);
        if isempty(h), target=pool(1); else, target=pool(h); end
    end
end

function v = canon_variant(v)
    v = string(v); L = lower(strtrim(v));
    v = strings(size(L));
    for i=1:numel(L)
        s = L(i);
        if any(strcmp(s,["original","ancestral","wuhan","wildtype","wild-type","wild type"])), v(i)="Alpha";
        elseif s=="alpha", v(i)="Alpha";
        elseif s=="delta", v(i)="Delta";
        elseif s=="omicron", v(i)="Omicron";
        else
            vv = upper(extractBefore(L(i),1))+extractAfter(L(i),1);
            if strlength(vv)==0, vv="Other"; end
            v(i)=vv;
        end
    end
end

function reg = region_series_of(T, colRegion)
    if strlength(colRegion)==0
        reg = repmat("Other", height(T),1);
    else
        reg = string(T.(colRegion));
    end
end

function R = province_region_mode(T, prov_levels, region_series)
    % 根据 Features 表里每条记录的 region，计算每省众数 region
    vn = string(T.Properties.VariableNames);
    hit = find(strcmpi(vn,"province") | strcmpi(vn,"Province"), 1);
    if isempty(hit), error('Features 表中未找到 province 列。'); end
    prov_col = vn(hit);
    R = strings(numel(prov_levels),1);
    for i = 1:numel(prov_levels)
        mask = string(T.(prov_col)) == prov_levels(i);
        regs = string(region_series(mask));
        if isempty(regs)
            R(i) = "Other";
        else
            u = unique(regs); cnt = zeros(size(u));
            for j = 1:numel(u), cnt(j)=sum(regs==u(j)); end
            [~, m] = max(cnt); R(i) = u(m);
        end
    end
end

function [mu_p, sd_p, ciP, mu_v, sd_v, var_levels] = read_ranef(root, target, prov_levels)
    % 省随机效应
    rpPath = fullfile(root,'M0_main',char(target),'ranef_province.csv');
    if ~isfile(rpPath), error('缺少文件：%s', rpPath); end
    rp = readtable(rpPath);
    [okP, iP] = ismember(prov_levels, string(rp.province));
    if ~all(okP)
        warning('部分省在 ranef_province.csv 中缺失：%s', strjoin(cellstr(prov_levels(~okP)),', '));
    end
    P = numel(prov_levels);
    mu_p = nan(P,1); sd_p = nan(P,1); ciP = nan(P,2);
    for i=1:P
        if iP(i)>0
            mu_p(i)=rp.mean(iP(i));
            sd_p(i)=rp.sd(iP(i));
            if ~(isfinite(sd_p(i))) && all(isfinite([rp.ci2_5(iP(i)) rp.ci97_5(iP(i))]))
                sd_p(i)=(rp.ci97_5(iP(i))-rp.ci2_5(iP(i)))/(2*1.96);
            end
            ciP(i,:)=[rp.ci2_5(iP(i)) rp.ci97_5(iP(i))];
        end
    end

    % 变异株随机效应（把 Original 合并为 Alpha）
    rvPath = fullfile(root,'M0_main',char(target),'ranef_variant.csv');
    if ~isfile(rvPath), error('缺少文件：%s', rvPath); end
    rv = readtable(rvPath);
    rv.variant = canon_variant(string(rv.variant));
    rv.variant = replace(rv.variant,"Original","Alpha");
    var_levels = unique(rv.variant,'stable');
    K = numel(var_levels);
    mu_v = nan(K,1); sd_v = nan(K,1);
    for k=1:K
        j = find(rv.variant==var_levels(k), 1);
        if ~isempty(j)
            mu_v(k) = rv.mean(j);
            sd_v(k) = rv.sd(j);
            if ~(isfinite(sd_v(k))) && all(isfinite([rv.ci2_5(j) rv.ci97_5(j)]))
                sd_v(k) = (rv.ci97_5(j)-rv.ci2_5(j))/(2*1.96);
            end
        end
    end
end

function W = compute_weights(T, prov_levels, var_levels, scheme, colCases, colDuration)
    P = numel(prov_levels); K = numel(var_levels);
    W = zeros(P,K);
    prov_ser = string(T.(pick(T,["province","Province"])));
    var_ser  = canon_variant(string(T.(pick(T,["variant","Variant"]))));
    var_ser  = replace(var_ser,"Original","Alpha");

    for i=1:P
        mP = prov_ser==prov_levels(i);
        for k=1:K
            mV = var_ser==var_levels(k);
            idx = find(mP & mV);
            if isempty(idx), val=0;
            else
                switch lower(scheme)
                    case 'cases'
                        if strlength(colCases)>0
                            x = double(T.(colCases)(idx)); x(~isfinite(x))=0; val = nansum(x);
                        else, val = numel(idx); end
                    case 'duration'
                        if strlength(colDuration)>0
                            x = double(T.(colDuration)(idx)); x(~isfinite(x))=0; val = nansum(x);
                        else, val = numel(idx); end
                    otherwise
                        val = numel(idx); % 波数
                end
            end
            W(i,k) = val;
        end
        s = sum(W(i,:)); if s>0, W(i,:) = W(i,:)/s; end
    end
end

function [blocks, order_idx, reg_names] = block_order_by_region(region_of_prov, mu_prov)
    [reg_names, ~, rid] = unique(region_of_prov,'stable');
    order_idx = []; blocks = [];
    for r=1:numel(reg_names)
        idx = find(rid==r);
        [~,ord] = sort(mu_prov(idx),'descend');
        order_idx = [order_idx; idx(ord)]; %#ok<AGROW>
        blocks    = [blocks;    r*ones(numel(idx),1)]; %#ok<AGROW>
    end
end

function r = nice_rmax(x)
    if isempty(x)||~isfinite(x), r=1; return; end
    a = max(0.2, x); pow = floor(log10(a)); base = a/(10.^pow);
    if base<=1, m=1; elseif base<=2, m=2; elseif base<=5, m=5; else, m=10; end
    r = m*(10.^pow);
end

function draw_circle(ax, r, col)
    th = linspace(0,2*pi,360); [x,y] = pol2cart(th, r);
    plot(ax, x, y, '-', 'Color', col, 'LineWidth',1.2, 'HandleVisibility','off');
end

function draw_wedge(ax, ang1, ang2, r1, r2, col, fa)
    if r2<r1, [r1,r2]=deal(r2,r1); end
    th = linspace(ang1, ang2, 40);
    [x1,y1] = pol2cart(th, r1);
    [x2,y2] = pol2cart(fliplr(th), r2);
    patch(ax, [x1 x2], [y1 y2], col, 'FaceAlpha', fa, 'EdgeColor','none', 'LineWidth',0.5);
end

function draw_mark(ax, ang, r, col)
    [x,y] = pol2cart(ang, r);
    plot(ax, x, y, 'o', 'MarkerSize', 5, 'MarkerFaceColor', col, ...
        'MarkerEdgeColor',[0.1 0.1 0.1], 'LineWidth',0.3);
end

function draw_separators(ax, theta, dtheta, blocks, rmax_plot)
    for k=2:numel(blocks)
        if blocks(k)~=blocks(k-1)
            ang = theta(k)-dtheta/2;
            [x1,y1] = pol2cart(ang, 0);
            [x2,y2] = pol2cart(ang, 1.12*rmax_plot);
            plot(ax, [x1 x2], [y1 y2], '-', 'Color',[0.82 0.82 0.82], 'LineWidth',0.8);
        end
    end
end

function col = region_color(Cgrp, rid, n_regions)
    if isempty(Cgrp)||size(Cgrp,2)~=3, Cgrp=lines(max(8,n_regions)); end
    K = size(Cgrp,1);
    if n_regions<=K
        idxs = round(linspace(1, K, n_regions)); idx = idxs(max(1, min(n_regions, rid)));
    else
        idx = mod(rid-1, K) + 1;
    end
    col = Cgrp(idx,:);
end

function v = safe_palette(f,fallback)
    try, v=f(); catch, v=fallback; end
end

function fn = pick_font(cands)
    try
        avail = listfonts;
        for ii = 1:numel(cands)
            if any(strcmpi(avail, cands{ii})), fn = cands{ii}; return; end
        end
        fn = avail{1};
    catch
        fn = 'Helvetica';  % 兜底
    end
end

end
