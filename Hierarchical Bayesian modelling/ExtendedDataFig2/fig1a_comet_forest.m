function fig = fig1a_comet_forest(results_root, varargin)
% FIG1A_COMET_FOREST — 固定效应“彗星森林图”（多目标并列）
% 显示每个目标（列）的 β 后验均值（彗星核）与 95%CI（渐变尾），并按变量组着色。
% 星标：基于正态近似的“后验同号概率” P_same = Φ(|mean/sd|)。
%
% 用法示例：
%   fig1a_comet_forest('results_hbm', ...
%       'palette_key','mainA8_deep', 'share_xlim',true, ...
%       'head_size',150, 'tail_width',7, 'nseg',9, 'dpi',300);

% ---------------- 参数 ----------------
p = inputParser;
addParameter(p,'max_rows',28,@(x)isnumeric(x)&&x>=5);
addParameter(p,'dpi',300,@(x)isnumeric(x)&&x>=72);

% 版式 & 视觉
addParameter(p,'font_tick',10,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'font_label',12,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'font_title',14,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'font_star',12,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'head_size',140,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'tail_width',6,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'nseg',8,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'tail_style','patch',@(s) any(strcmpi(s,{'taper','patch','flat'}))); % 推荐 patch
addParameter(p,'tail_min_ratio',0.25,@(x)isnumeric(x)&&x>0&&x<=1);   % r_min
addParameter(p,'tail_lambda',3.0,@(x)isnumeric(x)&&isscalar(x)&&x>0);% 形状λ，越大收得越快
addParameter(p,'tail_min_px',2.5,@(x)isnumeric(x)&&isscalar(x)&&x>0);% 尾部像素级最小厚度
addParameter(p,'tail_light_range',[0.15 0.70],@(v)isnumeric(v)&&numel(v)==2); % 渐亮范围
addParameter(p,'tail_light_lambda',2.2,@(x)isnumeric(x)&&isscalar(x)&&x>0);   % 颜色指数λ
addParameter(p,'tail_thick_data',0.30,@(x)isnumeric(x)&&isscalar(x)&&x>0);   % 基部厚度(数据单位)


% 配色与坐标口径
addParameter(p,'palette_key','groups8_deep',@(s)ischar(s)||isstring(s));
addParameter(p,'share_xlim',true,@islogical);
addParameter(p,'star_levels',[0.95 0.99 0.999],@(v)isnumeric(v)&&isvector(v));

parse(p,varargin{:});
args     = p.Results;
tick_fs  = args.font_tick; lab_fs = args.font_label; tit_fs = args.font_title;
star_fs  = args.font_star; HEAD_SZ= args.head_size;  TAIL_W = args.tail_width;
NSEG     = args.nseg;

if nargin<1||isempty(results_root), results_root='results_hbm'; end
S = load(fullfile(results_root,'cache','hbm_cache.mat')); cache = S.cache;
targets = string(fieldnames(cache.Y));
fe_names_full = string(cache.fe_names_orig);

% ---------------- 读各目标 summary.csv 并解析 beta 行 ----------------
Tbetas = cell(numel(targets),1);
for i=1:numel(targets)
    yname = targets(i);
    fsum = fullfile(results_root,'M0_main',yname,'summary.csv');
    if ~isfile(fsum), error('缺少 %s', fsum); end
    Ts = readtable(fsum,'TextType','string');
    Tbetas{i} = parse_beta_rows(Ts, fe_names_full);   % param, mean, sd, ci2_5, ci97_5
end

% ---------------- 统一变量集合 & 选 Top rows ----------------
param_cells = cellfun(@(Tb) string(Tb.param), Tbetas, 'UniformOutput', false);
all_names   = unique(vertcat(param_cells{:}), 'stable');

% 以“跨目标 |β| 的平均值”排序，取前 max_rows
metric = zeros(numel(all_names),1);
for j=1:numel(all_names)
    nm = all_names(j); acc=[];
    for i=1:numel(targets)
        Tb = Tbetas{i};
        hit = find(Tb.param==nm,1);
        if ~isempty(hit), acc(end+1)=abs(Tb.mean(hit)); end %#ok<AGROW>
    end
    metric(j) = mean(acc,'omitnan');
end
[~,ord]  = sort(metric,'descend');
keep     = ord(1:min(args.max_rows,numel(ord)));
row_names= all_names(keep);
nrows    = numel(row_names);
Y        = 1:nrows;  % 递增向量

% ---------------- 颜色映射（按变量组） ----------------
[group_of, cmap_of] = group_palette_with_other(args.palette_key);

% ---------------- 共享X轴（可选） ----------------
if args.share_xlim
    all_lo = []; all_hi = [];
    for i=1:numel(Tbetas)
        all_lo = [all_lo; Tbetas{i}.ci2_5]; %#ok<AGROW>
        all_hi = [all_hi; Tbetas{i}.ci97_5]; %#ok<AGROW>
    end
    xmax = max(abs([all_lo; all_hi]), [], 'omitnan');
    xlim_shared = [-1 1]*(xmax*1.05);  % 多 5% 余量
else
    xlim_shared = [];
end

% ---------------- 画图 ----------------
fig = figure('Name','Fig1a Comet Forest','Color','w','Position',[50 50 1200 800]);
t = tiledlayout(1, numel(targets),'Padding','compact','TileSpacing','compact');

for i = 1:numel(targets)
    nexttile; hold on;
    Tb = Tbetas{i};
    for r = 1:nrows
        nm  = row_names(r);
        hit = find(Tb.param==nm, 1);
        if isempty(hit), continue; end
        mu = Tb.mean(hit); lo = Tb.ci2_5(hit); hi = Tb.ci97_5(hit);
        sd = max(Tb.sd(hit), eps);

        gname = group_of(nm);
        col   = cmap_of(gname);

        draw_comet(mu, lo, hi, Y(r), col, HEAD_SZ, TAIL_W, NSEG, args);


        % 星标（后验同号概率）
        z = abs(mu/max(sd,eps));
        p_same = normcdf(z,0,1);
        levels = sort(args.star_levels(:)');
        nstars = sum(p_same >= levels);
        if nstars > 0
            text(mu, Y(r), [' ', repmat('*',1,nstars)], ...
                'Color', col*0.6, 'FontWeight','bold', 'FontSize', star_fs, ...
                'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
                'Interpreter','none');
        end
    end

    if ~isempty(xlim_shared), xlim(xlim_shared); end
    xline(0,'--','Color',[0.85 0.85 0.85],'LineWidth',0.6);

    % 轴与标签
    if i==1
        set(gca, 'YDir','reverse', 'YTick', Y, ...
            'YTickLabel', truncate_labels(row_names,30), ...
            'TickLabelInterpreter','none');
    else
        set(gca, 'YDir','reverse', 'YTick', Y, 'YTickLabel', repmat("", size(Y)), ...
            'TickLabelInterpreter','none');
    end
    title(string(targets(i)),'FontSize',tit_fs,'FontWeight','bold');
    xlabel('Effect size','FontSize',lab_fs);
    ylim([0 nrows+1]); box on; grid on;

    ax = gca;
    set(ax, 'FontName','Helvetica', 'FontSize', tick_fs, 'LineWidth', 0.6);
end

% hta = annotation(gcf,'textbox',[0.06 0.96 0.88 0.035], ... % [x y w h] 归一化
%     'String','Fixed effects across targets — Comet Forest', ...
%     'HorizontalAlignment','center','VerticalAlignment','middle', ...
%     'EdgeColor','none','FontSize',tit_fs+1,'FontWeight','bold', ...
%     'Interpreter','none','PickableParts','all','HitTest','on');

ht = sgtitle('Fixed effects across targets — Comet Forest');
set(ht,'FontSize',tit_fs+1,'FontWeight','bold','Interpreter','none');

% 兼容性：有啥属性就设啥（不同版本类名/属性不一致）
try
    if isprop(ht,'HitTest'),        set(ht,'HitTest','on');        end
    if isprop(ht,'PickableParts'),  set(ht,'PickableParts','all'); end
    if isprop(ht,'Pickable'),       set(ht,'Pickable','on');       end
catch
    % 忽略：老版本可能不支持这些属性，标题依然能编辑/导出
end



fig_make_editable(fig);
% 让文字/坐标轴可点，其它图元不抢点击
set(findall(fig,'Type','text'), 'PickableParts','all','HitTest','on');
set(findall(fig,'Type','axes'), 'Layer','top');
set(findall(fig,'Type','patch'),'HitTest','off','PickableParts','none');
set(findall(fig,'Type','line'), 'HitTest','off','PickableParts','none');
set(findall(fig,'Type','scatter'),'PickableParts','all','HitTest','on');


% ---------------- 导出 ----------------
outdir = fullfile(results_root,'M0_main','_figs');
if ~isfolder(outdir), mkdir(outdir); end
if args.dpi>0
    print(fig, fullfile(outdir,'Fig1a_comet_forest.png'), '-dpng', sprintf('-r%d',args.dpi));
end

end

% =========================================================================
% Helpers
% =========================================================================
function T = parse_beta_rows(Tsum, fe_names_orig)
% 把 summary.csv 中的 beta:* 与 beta_k 统一映射为“可读名字”
is_named = startsWith(Tsum.param,"beta:");
is_idxed = ~is_named & startsWith(Tsum.param,"beta_");
A = table(); B = table();

if any(is_named)
    Tn = Tsum(is_named,:);
    A = table(extractAfter(Tn.param,"beta:"), Tn.mean, Tn.sd, Tn.ci2_5, Tn.ci97_5, ...
        'VariableNames',{'param','mean','sd','ci2_5','ci97_5'});
end
if any(is_idxed)
    Ti = Tsum(is_idxed,:);
    k  = double(extractAfter(Ti.param,"beta_"));
    nm = strings(height(Ti),1);
    for r=1:height(Ti)
        if ~isnan(k(r)) && k(r)>=1 && k(r)<=numel(fe_names_orig)
            nm(r) = string(fe_names_orig{k(r)});
        else
            nm(r) = string(Ti.param(r));
        end
    end
    B = table(nm, Ti.mean, Ti.sd, Ti.ci2_5, Ti.ci97_5, ...
        'VariableNames',{'param','mean','sd','ci2_5','ci97_5'});
end
T = [A;B];
[~,ia] = unique(T.param,'stable'); T = T(ia,:);  % 去重复，保留首次
end

function draw_comet(mu, lo, hi, y, col, markerSize, tailWidth, nseg, args)
% 指数型渐细彗尾 + 像素级最小厚度保护
if nargin < 6 || isempty(markerSize), markerSize = 140; end
if nargin < 7 || isempty(tailWidth),  tailWidth  = 6;   end
if nargin < 8 || isempty(nseg),       nseg = 10;        end
if nargin < 9, args = struct(); end

style    = lower(args.tail_style);              % 'patch' | 'taper' | 'flat'
rmin     = args.tail_min_ratio;                 % 尾端/基部 相对厚度
lambda   = args.tail_lambda;                    % 几何收束 λ
min_px   = args.tail_min_px;                    % 像素下限
w0w1     = args.tail_light_range(:)';           % 渐亮范围 [w0 w1]
lam_c    = args.tail_light_lambda;              % 颜色指数 λ
baseH    = args.tail_thick_data;                % 基部厚度(数据单位)
bg       = [1 1 1];

% --- 工具：颜色向白指数混合 (w(t)指数增长) ---
mix_white = @(c,t) c*(1-(w0w1(1)+(w0w1(2)-w0w1(1))*(1-exp(-lam_c*t)))) ...
                     + bg*(w0w1(1)+(w0w1(2)-w0w1(1))*(1-exp(-lam_c*t)));

% --- 计算像素→数据单位的换算，保障最小厚度 ---
ax = gca;
axpix = getpixelposition(ax,true);       % [x y w h] in pixels
ylimv = ylim(ax);
dy_per_px = diff(ylimv)/max(axpix(4),1); % 每像素对应的y数据单位
minH_data = max(min_px*dy_per_px, rmin*baseH*0.9);  % 两个下限取较大者

% --- 指数厚度函数（t∈[0,1]，近核t小、远端t大） ---
hfun = @(t) max( baseH*( rmin + (1-rmin)*exp(-lambda*t) ), minH_data );

switch style
    case 'patch'
        % 左尾（mu ← lo）
        if lo < mu
            xs = linspace(mu, lo, nseg+1);
            for k=1:nseg
                t1 = (k-1)/nseg; t2 = k/nseg;
                h1 = hfun(t1);   h2 = hfun(t2);
                colk = mix_white(col, t2);
                x1 = xs(k); x2 = xs(k+1);
                X  = [x1 x2 x2 x1];
                Y  = [y-h1/2, y-h2/2, y+h2/2, y+h1/2];
                patch(X, Y, colk, 'EdgeColor','none');
            end
        end
        % 右尾（mu → hi）
        if hi > mu
            xs = linspace(mu, hi, nseg+1);
            for k=1:nseg
                t1 = (k-1)/nseg; t2 = k/nseg;
                h1 = hfun(t1);   h2 = hfun(t2);
                colk = mix_white(col, t2);
                x1 = xs(k); x2 = xs(k+1);
                X  = [x1 x2 x2 x1];
                Y  = [y-h1/2, y-h2/2, y+h2/2, y+h1/2];
                patch(X, Y, colk, 'EdgeColor','none');
            end
        end

    case 'taper'
        % 线段指数收细（轻量）：LineWidth 用 hfun(t) 换算成“等效粗度”
        scaleLW = 1/0.75; % 调整：把数据厚度映射到合适的LW视觉量级
        if lo < mu
            xs = linspace(mu, lo, nseg+1);
            for k=1:nseg
                t  = k/nseg;
                lw = max( tailWidth * ( hfun(t)/baseH )*scaleLW, 0.5 );
                colk = mix_white(col, t);
                plot([xs(k+1) xs(k)], [y y], 'Color', colk, 'LineWidth', lw);
            end
        end
        if hi > mu
            xs = linspace(mu, hi, nseg+1);
            for k=1:nseg
                t  = k/nseg;
                lw = max( tailWidth * ( hfun(t)/baseH )*scaleLW, 0.5 );
                colk = mix_white(col, t);
                plot([xs(k) xs(k+1)], [y y], 'Color', colk, 'LineWidth', lw);
            end
        end

    otherwise % 'flat' 维持旧版
        if lo < mu
            xs = linspace(mu, lo, nseg+1);
            for k=1:nseg
                t = k/nseg; colk = mix_white(col, t);
                plot([xs(k+1) xs(k)], [y y], 'Color', colk, 'LineWidth', tailWidth);
            end
        end
        if hi > mu
            xs = linspace(mu, hi, nseg+1);
            for k=1:nseg
                t = k/nseg; colk = mix_white(col, t);
                plot([xs(k) xs(k+1)], [y y], 'Color', colk, 'LineWidth', tailWidth);
            end
        end
end

% —— 彗星核 —— %
edge = max(col*0.6, 0);
scatter(mu, y, markerSize, 'MarkerFaceColor', col, 'MarkerEdgeColor', edge, 'LineWidth', 0.8);
end



function s = truncate_labels(a, maxch)
s = strings(size(a));
for i=1:numel(a)
    t = char(a(i)); if numel(t)>maxch, t = [t(1:maxch-1) '…']; end
    s(i) = string(t);
end
end

function [group_of, cmap_of] = group_palette_with_other(palette_key)
% 8 组关键词映射（与论文蓝图一致）
grp = struct( ...
    'Population', ["age" "student" "education"], ...
    'Healthcare', ["hospital" "beds"], ...
    'Economy',    ["gdp" "income"], ...
    'Urbanicity', ["urban" "density" "tertiary"], ...
    'Mobility',   ["passenger" "transit"], ...
    'Tourism',    ["inbound" "overnight" "visitors"], ...
    'Region',     ["region"], ...
    'Climate',    ["climate"] ...
    );

% 从色库拿 8 色；失败则回退到 groups8_deep
try
    C8 = smf_colorlib(palette_key,'K',8);
catch
    C8 = smf_colorlib('groups8_deep','K',8);
end
order = ["Population","Healthcare","Economy","Urbanicity","Mobility","Tourism","Region","Climate"];
cmap  = containers.Map(order, num2cell(C8,2));
cmap('Other') = [0.55 0.55 0.55];

keys = fieldnames(grp);
group_of = @(nm) local_map_group(lower(string(nm)), grp, keys);
cmap_of  = @(g) cmap_get(g);

    function g = local_map_group(nm, grp, keys)
        for i=1:numel(keys)
            K = keys{i}; toks = grp.(K);
            for j=1:numel(toks)
                if contains(nm, toks(j)), g = string(K); return; end
            end
        end
        g = "Other";
    end
    function col = cmap_get(g)
        k = char(g);
        if isKey(cmap, k), col = cmap(k); else, col = cmap('Other'); end
    end
end

function fig_make_editable(fig)
if nargin<1 || ~ishandle(fig), fig = gcf; end
set(findall(fig,'Type','text'),    'PickableParts','all','HitTest','on');
set(findall(fig,'Type','axes'),    'Layer','top');
set(findall(fig,'Type','scatter'), 'PickableParts','all','HitTest','on');
set(findall(fig,'Type','patch'),   'HitTest','off','PickableParts','none');
set(findall(fig,'Type','line'),    'HitTest','off','PickableParts','none');
% sgtitle/annotation 也设置一下
h = findall(fig,'Type','text','-and','Tag','Sgt'); % 某些版本sgtitle带Tag
set(h,'PickableParts','all','HitTest','on');
end
