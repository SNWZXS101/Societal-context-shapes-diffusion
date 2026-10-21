%% ================= FIG.2c1 — Variant Composition (Bar + Donuts, RUN) =================
% 上：各变异株(A/D/O) 内 C1–C4 并列条 + 95% CI + 省级占比散点（半透明）
% 下：每个变异株一张“环形饼图”（donut），展示该毒株内 C1–C4 的占比（直观对比）
% 统计：全局 χ² (variant×cluster) 与 Cramér’s V；每个变异株内对“均匀分布”的χ² p
% 字号统一 12；标题为 annotation，可在图窗里直接编辑
clear; clc; close all;

%% ------------------- 用户配置 -------------------
results_mat   = 'smf_fpca_cluster_results.mat';
features_mat  = 'features.mat';
outdir        = 'figs_main';  if ~exist(outdir,'dir'), mkdir(outdir); end
B_boot        = 2000;  rng(1);           % bootstrap 次数
USE_WEIGHT    = false;                   % 散点大小是否用省内样本数（条/CI仍按计数）
VARIANT_ORDER = ["Alpha","Delta","Omicron"];   % 显示顺序（缺谁跳谁）
DOT_ALPHA     = 0.55;                    % 散点透明度
JITTER_W      = 0.02;                    % 上图散点的最大抖动
FS            = 15;                      % 统一字号
% —— p 括号位置控制（相对 y 轴范围的百分比；也可以给每组单独值）——
PBRACKET_Y_FRAC = 0.96;   % 括号的基准高度 = ymin + 0.96*(ymax-ymin)
PBRACKET_STEP   = 0.05;   % 多个括号时的上下间距（乘以 y-range）
PBRACKET_TEXT_DY= 0.015;  % p 值文字在括线上方再抬一点（乘以 y-range）


% —— 上图条形的“组间距/组内距”控制（重点：用 INNER_PACK）——
GROUP_GAP   = 1.0;   % 组中心间隔的单位刻度（保持 1 即可）
INNER_PACK  = 0.50;  % 组内条整体占比（越小→组间留白越大；建议 0.55~0.85）
F_BAR       = 0.9;  % 单条在槽位中的占比（越小→条更瘦；建议 0.6~0.85）

% —— 下图 donut 的外观 ——
DONUT_INNER_FRAC = 0.55;   % 内孔相对半径（0.45~0.65 合适）
DONUT_PAD        = 0.04;   % 片间小缝（弧度），避免色块黏连

%% ------------------- 全局字体 -------------------
set(groot,'defaultAxesFontSize',FS);
set(groot,'defaultTextFontSize',FS);
set(groot,'defaultLegendFontSize',FS);

%% ------------------- 读取/规整 -------------------
S = load(results_mat);
[labels, K] = get_labels_from_results(S);         % K=4
[variant, province, w] = get_meta_variant_prov_weight(S, features_mat, numel(labels));
if isempty(w), w = ones(numel(labels),1); end

uv  = unique(variant,'stable');
ord = VARIANT_ORDER(ismember(VARIANT_ORDER, uv));
if isempty(ord), ord = uv; end
G = numel(ord);

% 变异株×簇计数/占比
Tab_vc = zeros(G, K);
for i=1:G
    vi = (variant==ord(i));
    for k=1:K, Tab_vc(i,k) = sum(vi & labels==k); end
end
RowSum = sum(Tab_vc,2); RowSum(RowSum==0)=1;
Prop_vc = Tab_vc ./ RowSum;

% 全局 χ² 与 V
[p_glb, chi2_glb, df_glb] = chisq_independence(Tab_vc);
V_glb = cramers_v(chi2_glb, sum(Tab_vc,'all'), size(Tab_vc,1), size(Tab_vc,2));

% 组内对“均匀”的 χ² p
p_group = nan(G,1);
for i=1:G
    Oi = Tab_vc(i,:); ni = sum(Oi);
    if ni>0
        Ei = ones(1,K) * (ni/K);
        p_group(i) = chisq_pvalue(Oi, Ei);
    end
end

% 颜色
Cclus = smf_colorlib('root');         % C1–C4 色

%% ------------------- 省级占比（供上图散点） -------------------
ProvShare = cell(G, K);
ProvN     = cell(G, 1);
for i=1:G
    idx_v = (variant==ord(i));
    prov_i = unique(province(idx_v),'stable');
    Pi = numel(prov_i);
    ProvN{i} = zeros(Pi,1);
    for pidx=1:Pi
        idx_pv = idx_v & (province==prov_i(pidx));
        ProvN{i}(pidx) = sum(idx_pv);
        li = labels(idx_pv);
        if isempty(li), continue; end
        for k=1:K
            ProvShare{i,k}(end+1,1) = mean(li==k); %#ok<AGROW>
        end
    end
end

%% ------------------- 画图容器 -------------------
f = figure('Color','w','Position',[60 60 1150 820]);
TL = tiledlayout(f,2,1,'TileSpacing','compact','Padding','compact');

%% ===== 上面板：条形 + 95%CI + 省级散点（组间距由 INNER_PACK 控）=====
ax1 = nexttile(TL,1); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');

% 组中心（用单位间隔 1 的刻度，视觉组间留白由 INNER_PACK 决定）
xg = (1:G) * GROUP_GAP;                 % 实际上就是 1:G
inner_span = GROUP_GAP * INNER_PACK;    % 组内总宽度（相对组中心刻度）
offsets    = linspace(-inner_span/2, inner_span/2, K);
bar_abs_w  = (inner_span / K) * F_BAR;
bar_rel_w  = bar_abs_w / GROUP_GAP;     % bar() 需要相对宽度

% 逐簇绘制（而非 'grouped'）
Ybar = Prop_vc; b = gobjects(1,K);
for k = 1:K
    xk = xg + offsets(k);
    b(k) = bar(ax1, xk, Ybar(:,k), ...
        'BarWidth',  bar_rel_w, ...
        'FaceColor', Cclus(k,:), ...
        'EdgeColor', [0.2 0.2 0.2], ...
        'LineWidth', 0.5);
end
ylim(ax1,[0 1]); ylabel(ax1,'Proportion within variant');
set(ax1,'XTick',xg,'XTickLabel',cellstr(ord));

% 组分隔竖线（淡灰）
yl = ylim(ax1);
for i=1:G-1
    xsep = (xg(i)+xg(i+1))/2;
    plot(ax1,[xsep xsep], yl, 'Color',[0.88 0.88 0.88], 'LineStyle','-', 'LineWidth',0.8, 'HandleVisibility','off');
end

% 每根柱中心 x（CI/散点使用）
xCenters = repmat(xg(:),1,K) + repmat(offsets(:).',G,1);

% bootstrap 95% CI
for i=1:G
    idx_v = (variant==ord(i)); lab_v = labels(idx_v);
    if numel(lab_v)<=1, continue; end
    Pboot = zeros(B_boot, K);
    for bb=1:B_boot
        samp = lab_v(randi(numel(lab_v), numel(lab_v), 1));
        for k=1:K, Pboot(bb,k) = mean(samp==k); end
    end
    lo = prctile(Pboot, 2.5, 1); hi = prctile(Pboot, 97.5, 1);
    for k=1:K
        xc = xCenters(i,k);
        plot(ax1, [xc xc], [lo(k) hi(k)], '-', 'Color',[0 0 0], 'LineWidth',0.8, 'HandleVisibility','off');
        plot(ax1, [xc-0.05 xc+0.05], [lo(k) lo(k)], '-', 'Color',[0 0 0], 'LineWidth',0.8, 'HandleVisibility','off');
        plot(ax1, [xc-0.05 xc+0.05], [hi(k) hi(k)], '-', 'Color',[0 0 0], 'LineWidth',0.8, 'HandleVisibility','off');
    end
end

% 省级散点（抖动不超过“邻条中心距离的一半”）
for i=1:G
    for k=1:K
        if isempty(ProvShare{i,k}), continue; end
        xc = xCenters(i,k);
        neigh = sort(xCenters(i,:));
        step  = min(diff(neigh));           % 邻条中心距离
        jitter_max = 0.45 * step;           % 安全抖动上限
        xj = xc + (rand(size(ProvShare{i,k}))-0.5) * 2 * min(JITTER_W, jitter_max);
        scatter(ax1, xj, ProvShare{i,k}, 26, Cclus(k,:), 'filled', ...
                'MarkerFaceAlpha', DOT_ALPHA, 'MarkerEdgeAlpha', 0, ...
                'HandleVisibility','off');
    end
end

% 每组上方 p 值括号（括住该组4根柱的最左/最右中心）
yl = ylim(ax1); yr = yl(2)-yl(1);
HB = cell(G,1);
for i=1:G
    xl = min(xCenters(i,:)); xr = max(xCenters(i,:));
    yb = yl(1) + (PBRACKET_Y_FRAC - (i-1)*PBRACKET_STEP) * yr;
    pstr = sprintf('p %s', pretty_p(p_group(i)));
    HB{i} = draw_pbracket_anno(ax1, xl, xr, yb, pstr);  % ← 现在是annotation，可拖动
end
plotedit on;  % 开启图窗交互编辑（可选）


% 组下方 n
for i=1:G
    text(ax1, mean(xCenters(i,:)), -0.08, sprintf('n=%d', sum(variant==ord(i))), ...
        'HorizontalAlignment','center');
end

legend(ax1, arrayfun(@(k)sprintf('C%d',k),1:K,'UniformOutput',false), ...
       'Location','eastoutside');

%% ===== 下面板：每个变异株一张“环形饼图”展示 C1–C4 占比（直观）=====
ax_dummy = nexttile(TL,2); axis(ax_dummy,'off');  % 占位，不直接用
delete(ax_dummy);

% 手动创建三块小轴，水平排布
margL=0.08; margR=0.08; margT=0.06; margB=0.10; gap=0.04;
w = (1 - margL - margR - (G-1)*gap)/G; h = 0.34;  % 下排高度
y0 = margB;
for i=1:G
    x0 = margL + (i-1)*(w+gap);
    axD(i) = axes('Position',[x0 y0 w h]); %#ok<AGROW>
    hold(axD(i),'on'); axis(axD(i),'equal'); axis(axD(i),'off');

    % 数据与颜色
    p = Prop_vc(i,:);                      % 该毒株内 C1..C4 占比
    draw_donut(axD(i), p, Cclus, DONUT_INNER_FRAC, DONUT_PAD);

    % 中心数值（可选：写上最大份额及其簇号）
    [mx, ik] = max(p);
    text(0,0,sprintf('C%d\n%.0f%%',ik,100*mx), ...
        'HorizontalAlignment','center','FontWeight','bold');

    % 小标题：变异株 + 样本数
    title(axD(i), sprintf('%s (n=%d)', ord(i), sum(variant==ord(i))));
end

%% ------------------- 可编辑主标题（annotation） -------------------
mainTitle = sprintf('Variant composition — global \\chi^2 p=%.3g, V=%.2f, N=%d', ...
                     p_glb, V_glb, numel(labels));
ht = annotation(f,'textbox',[0.12 0.94 0.76 0.05], 'String', mainTitle, ...
     'EdgeColor','none','HorizontalAlignment','center', ...
     'FontWeight','bold','FontSize',FS,'Interpreter','tex','FitBoxToText','off');
plotedit on;

%% ------------------- 保存 -------------------
cc = call_print_counts(labels, 'Tag','Fig2c1', 'ExpectK',4);
cc_str = strjoin(string(cc),'-');
fname = fullfile(outdir, sprintf('fig2c1_variant_bar_donut_K%d_[%s]', K, cc_str));
print(f,[fname,'.png'],'-dpng','-r300'); savefig(f,[fname,'.fig']);
disp(['Saved: ', fname, '.png & .fig']);

%% ====================== 辅助函数们（本脚本内置） ======================
function [labels, K] = get_labels_from_results(S)
    flds = fieldnames(S);
    if numel(flds)==1 && isstruct(S.(flds{1})), S=S.(flds{1}); end
    assert(isfield(S,'cluster') && isfield(S.cluster,'labels'), 'results_mat 缺少 cluster.labels');
    labels = S.cluster.labels(:); K = max(labels);
end
function [variant, province, w] = get_meta_variant_prov_weight(S, features_mat, N)
    variant = strings(0,1); province = strings(0,1); w = [];
    M=[]; if isfield(S,'meta'), M=S.meta; elseif isfield(S,'S')&&isfield(S.S,'meta'), M=S.S.meta; end
    pick = @(T, keys) local_pick(T, keys);
    if ~isempty(M)
        variant  = tostr(pick(M,["variant","Variant","strain"]));
        province = tostr(pick(M,["province","Province","region","prov"]));
        w        = pick(M,["weight","w","omega","weights"]);
    end
    if (numel(variant)~=N || numel(province)~=N || isempty(w)) && exist(features_mat,'file')
        T = load(features_mat); Tab = pick_widest_table(T);
        if ~isempty(Tab)
            if any(strcmpi(Tab.Properties.VariableNames,'variant')),  variant  = tostr(Tab.variant); end
            if any(strcmpi(Tab.Properties.VariableNames,'province')), province = tostr(Tab.province); end
            if isempty(province) && any(strcmpi(Tab.Properties.VariableNames,'prov')), province = tostr(Tab.prov); end
            if isempty(w)
                if any(strcmpi(Tab.Properties.VariableNames,'w_tail_14')), w = Tab.w_tail_14;
                elseif any(strcmpi(Tab.Properties.VariableNames,'w_tail')),  w = Tab.w_tail;
                elseif any(strcmpi(Tab.Properties.VariableNames,'weight')),  w = Tab.weight; end
            end
        end
    end
    if isempty(variant),  variant  = repmat("NA", N,1); end
    if isempty(province), province = repmat("NA", N,1); end
    if isempty(w),        w = ones(N,1); else, w = w(:); end
    variant  = tostr(variant); province = tostr(province);
end
function val = local_pick(T, keys)
    val=[]; if isempty(T), return; end
    if istable(T)
        for i=1:numel(keys), k=char(keys(i));
            if any(strcmpi(T.Properties.VariableNames,k)), val=T.(k); return; end
        end
    elseif isstruct(T)
        f=fieldnames(T);
        for i=1:numel(keys), k=char(keys(i));
            hit=f(strcmpi(f,k)); if ~isempty(hit), val=T.(hit{1}); return; end
        end
    end
end
function Tab = pick_widest_table(T)
    Tab=[]; fn=fieldnames(T); bestW=-1; key='';
    for i=1:numel(fn)
        if istable(T.(fn{i})) && width(T.(fn{i}))>bestW
            bestW=width(T.(fn{i})); key=fn{i};
        end
    end
    if ~isempty(key), Tab=T.(key); end
end
function s = tostr(x)
    if isstring(x), s=x; elseif ischar(x), s=string(x);
    elseif iscellstr(x), s=string(x);
    elseif iscell(x), try s=string(x); catch, s=strings(numel(x),1); for i=1:numel(x), if ischar(x{i}), s(i)=string(x{i}); end, end, end
    else, s=string(x); end
end

% —— 统计 —— 
function [p,chi2_stat,df] = chisq_independence(T)
    RT = sum(T,2); CT = sum(T,1); N = sum(RT);
    E  = (RT * CT) / max(N, eps);
    valid = E > 0;
    chi2_stat = sum(((T(valid)-E(valid)).^2) ./ E(valid), 'all');
    df = (size(T,1)-1)*(size(T,2)-1);
    p  = 1 - gammainc(chi2_stat/2, df/2, 'lower');
end
function V = cramers_v(chi2_stat, N, r, c)
    V = sqrt(chi2_stat / (N * max(1, min(r-1, c-1))));
end
function p = chisq_pvalue(O, E)
    valid = E>0;
    chi2 = sum(((O(valid)-E(valid)).^2) ./ E(valid), 'all');
    df   = sum(valid)-1;
    p    = 1 - gammainc(chi2/2, df/2, 'lower');
end

% —— p 值括号
function draw_pbracket(ax, x0, x1, y, p)
    if isnan(p), return; end
    plot(ax, [x0 x1], [y y], 'k-', 'LineWidth',0.8, 'HandleVisibility','off');
    plot(ax, [x0 x0], [y-0.01 y], 'k-', 'LineWidth',0.8, 'HandleVisibility','off');
    plot(ax, [x1 x1], [y-0.01 y], 'k-', 'LineWidth',0.8, 'HandleVisibility','off');
    t = sprintf('p %s', pretty_p(p));
    text(ax, mean([x0 x1]), y+0.015, t, 'HorizontalAlignment','center', 'VerticalAlignment','bottom');
end

% —— 把轴数据坐标转成 figure 的归一化坐标（用于 annotation 定位）——
function [xf, yf] = axdata2fig(ax, x, y)
    % 轴坐标→归一化 figure 坐标
    fig = ancestor(ax,'figure');
    axpos = getpixelposition(ax, true);  % [x y w h] in pixels relative to figure
    figpos = getpixelposition(fig, true);
    xl = xlim(ax); yl = ylim(ax);
    xp = axpos(1)/figpos(3) + (x - xl(1)) / (xl(2)-xl(1)) * (axpos(3)/figpos(3));
    yp = axpos(2)/figpos(4) + (y - yl(1)) / (yl(2)-yl(1)) * (axpos(4)/figpos(4));
    xf = xp; yf = yp;
end

% —— 用 annotation 画“可拖动”的 p 括号（返回句柄便于后续微调/删除）——
function h = draw_pbracket_anno(ax, x0, x1, y, pstr)
    fig = ancestor(ax,'figure');
    % 竖线高度（括号两端的小竖线），和横线 y 相同
    [x0f, yf] = axdata2fig(ax, x0, y);
    [x1f, ~ ] = axdata2fig(ax, x1, y);
    y0f = yf - 0.02;  % 竖线向下 2% fig 高度（可以改）
    % 左竖线、横线、右竖线
    h.L = annotation(fig,'line', [x0f x0f], [y0f yf], 'Color','k','LineWidth',0.8);
    h.H = annotation(fig,'line', [x0f x1f], [yf  yf],  'Color','k','LineWidth',0.8);
    h.R = annotation(fig,'line', [x1f x1f], [y0f yf], 'Color','k','LineWidth',0.8);
    % 文本（可拖动）
    h.T = annotation(fig,'textbox', [mean([x0f x1f])-0.03, yf+0.005, 0.08, 0.03], ...
        'String', pstr, 'FitBoxToText','on', 'EdgeColor','none', ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom');
end


function s = pretty_p(p)
    if p < 1e-3, s = '< 0.001';
    else, s = sprintf('= %.3f', p);
    end
end

% —— Donut 绘制（无依赖）——
function draw_donut(ax, propRow, C, innerFrac, pad)
    p = propRow(:)'; p = p / max(sum(p),eps);
    ang = [0, cumsum(p)*2*pi];
    R  = 1; r = innerFrac;
    for i=1:numel(p)
        a0 = ang(i)   + pad/2;
        a1 = ang(i+1) - pad/2;
        if a1<=a0, continue; end
        tt = linspace(a0,a1,60);
        x1 = R*cos(tt); y1 = R*sin(tt);
        x2 = r*cos(fliplr(tt)); y2 = r*sin(fliplr(tt));
        patch(ax, [x1 x2], [y1 y2], C(i,:), 'EdgeColor','none');
    end
    % 外圈淡灰描边（更精致）
    th = linspace(0,2*pi,200);
    plot(ax, R*cos(th), R*sin(th), 'Color',[0.85 0.85 0.85], 'LineWidth',0.8);
    plot(ax, r*cos(th), r*sin(th), 'Color',[1 1 1]*0.98, 'LineWidth',0.8);
end

% —— 配色兜底
% function C = smf_colorlib_safe(which)
%     if exist('smf_colorlib','file')==2
%         C = smf_colorlib(which);
%         return;
%     end
%     switch lower(which)
%         case {'classic4','clusters4','clusters4_deep'}
%             C = [0.12 0.47 0.71; 0.17 0.63 0.17; 0.58 0.40 0.74; 0.84 0.15 0.16];
%         case {'variants3','variants3_deep'}
%             C = [0.12 0.47 0.71; 1.00 0.50 0.05; 0.17 0.63 0.17];
%         otherwise
%             C = lines(4);
%     end
% end

% —— 计数打印（若全局有 smf_print_cluster_counts.m 就用它，否则用本地兜底）
function counts = call_print_counts(labels, varargin)
    if exist('smf_print_cluster_counts','file')==2
        counts = smf_print_cluster_counts(labels, varargin{:});
    else
        counts = smf_print_cluster_counts_local(labels, varargin{:});
    end
end
function counts = smf_print_cluster_counts_local(labels, varargin)
    p = inputParser;
    p.addParameter('Tag','',@(x)ischar(x)||isstring(x));
    p.addParameter('ExpectK',[],@(x)isempty(x)||(isscalar(x)&&x==round(x)));
    p.parse(varargin{:});
    tag = char(p.Results.Tag); Kexp = p.Results.ExpectK;

    labels = labels(:);
    if any(~isfinite(labels)), warning('labels 中存在 NaN/Inf，将被忽略计数。'); end
    L = labels(isfinite(labels));
    if isempty(L), fprintf('[%s] labels 为空。\n', tag); counts=[]; return; end
    K = max(L);
    counts = accumarray(max(1,min(K, round(L))), 1, [K 1])';
    N = numel(L);
    if ~isempty(tag), tinfo=[' [' tag ']']; else, tinfo=''; end
    fprintf('%s Cluster counts (K=%d, N=%d): %s\n', tinfo, K, N, mat2str(counts));
    if ~isempty(Kexp) && K ~= Kexp, warning('%s 实际 K=%d 与期望 K=%d 不一致。', tinfo, K, Kexp); end
    if sum(counts) ~= N, warning('%s 计数之和(%d) 与 N(%d) 不一致。', tinfo, sum(counts), N); end
    for k = 1:K, fprintf('  C%d: %d\n', k, counts(k)); end
end
