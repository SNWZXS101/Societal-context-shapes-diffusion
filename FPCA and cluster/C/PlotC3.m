%% ===== FIG.2c3 — Province → Cluster chord（RUN，增强版） =====
% - 左：四个簇 C1–C4；右：Top-N 省份（其余为 "Other"）
% - 宽弧带节点 + 贝塞尔丝带；左右弧段比例可配（默认 1/5 vs 3/5），双侧留白
% - 自动缩放 + 自动扩坐标以包含所有标签；省名自动截断；对象可选中编辑
clear; clc; close all;

%% 配置
results_mat     = 'smf_fpca_cluster_results.mat';
features_mat    = 'features.mat';
outdir          = 'figs_main';
TOPN_PROV       = 30;          % 右侧展示的省份数（其余合并 "Other"）
USE_WEIGHT      = false;       % false=计数；true=权重和（w_tail_14）

% —— 布局比例与间隙 ——
LEFT_FRAC       = 1/5;         % 原型簇占整圆比例
RIGHT_FRAC      = 3/5;         % 省份占整圆比例
INNER_GAP_FRAC  = 0.005;       % 组内相邻节点的细小间隙（占整圆比例）
ROT_DEG         = -315;        % 整图旋转角（度），正=逆时针，负=顺时针
% 其余 1 - LEFT_FRAC - RIGHT_FRAC 作为“组间留白”，自动均分到左右两端

% —— 尺寸控制 ——
SCALE           = 1;        % 整体半径缩放（0~1），让图更紧凑
R_OUTER         = 1.00 * SCALE;    % 外半径（节点弧外沿 / 标签参考）
R_INNER         = 0.9 * SCALE;    % 丝带起止半径（略小于节点弧内沿）
ARC_THICK       = 0.06;            % 节点弧带厚度（0.06~0.10 建议）
R_NODE_OUT      = R_OUTER;         % 节点弧外沿
R_NODE_IN       = R_NODE_OUT - ARC_THICK;    % 节点弧内沿
R_CTRL          = 0.30 * SCALE;    % 丝带贝塞尔控制半径（越小越弯）
ALPHA_RIBBON    = 0.85;            % 丝带透明度

% —— 标签样式 ——
LABEL_OFFSET    = 0.06;        % 标签相对外半径的外扩距离
MAX_LABEL_LEN   = 12;          % 省名最长字符数，超出截断加省略号
FONTSZ_LABEL    = 11;

SEED            = 1; rng(SEED);
if ~exist(outdir,'dir'), mkdir(outdir); end

%% 读取数据
S = load(results_mat);
[labels, K] = get_labels_from_results(S);
[variant, province, w] = get_meta_variant_prov_weight(S, features_mat, numel(labels));
assert(K==4, '当前脚本默认 K=4（如需改请相应调整簇颜色与标注）。');

% 省份 Top-N
prov_all = unique(province,'stable');
counts_by_prov = arrayfun(@(p) sum(province==p), prov_all);
[~, ord] = sort(counts_by_prov,'descend');
prov_top = prov_all(ord(1:min(TOPN_PROV, numel(prov_all))));
if numel(prov_top) < numel(prov_all), prov_top(end+1) = "Other"; end

% 省份×簇 流量矩阵 F (P×K)
P = numel(prov_top);
F = zeros(P, K);
for i=1:P
    if prov_top(i)=="Other"
        idxp = ~ismember(province, prov_top(1:end-1));
    else
        idxp = (province==prov_top(i));
    end
    for k=1:K
        F(i,k) = USE_WEIGHT * sum(w(idxp & labels==k), 'omitnan') + ...
                (~USE_WEIGHT) * sum(idxp & labels==k);
    end
end
if sum(F,'all')==0, warning('无可视化流量，退出'); return; end


cc = smf_print_cluster_counts(labels, 'Tag','Fig2c3', 'ExpectK',4);
cc_str = strjoin(string(cc), '-');


%% 配色
% 簇色建议用你论文主配色：clusters4_deep。如果你偏爱 root 也可保留。
% Cclus = smf_colorlib('clusters4_deep');
Cclus = smf_colorlib('root');           % ← 你提供的版本里用的是 'root'
Cprov = repmat([0.75 0.75 0.75], P, 1); % 右半：省份统一灰（也可换为分区配色）

%% 节点弧布局（整圆 + 双侧留白 + 可旋转）
left_total  = sum(F,1);         % 1×K
right_total = sum(F,2).';       % 1×P
assert(abs(sum(left_total) - sum(right_total)) < 1e-8, '左右总流量不等。');

gap_item = 2*pi*INNER_GAP_FRAC;
ANG0     = pi/2;                        % 从 12 点方向开始
SPAN_L   = 2*pi*LEFT_FRAC;
SPAN_R   = 2*pi*RIGHT_FRAC;
FREE_FRAC = max(0, 1 - LEFT_FRAC - RIGHT_FRAC);
gap_group_each = pi * FREE_FRAC;        % 左右各一半留白

% 左组（簇）
startL = ANG0 + gap_group_each;
theta_left  = alloc_arcs(left_total,  startL + gap_item/2,  SPAN_L - gap_item, gap_item);
% 右组（省）
startR = startL + SPAN_L + gap_group_each;
theta_right = alloc_arcs(right_total, startR + gap_item/2,  SPAN_R - gap_item, gap_item);

% 整体旋转
ROT_ANG = deg2rad(ROT_DEG);
theta_left  = rot_theta(theta_left,  ROT_ANG);
theta_right = rot_theta(theta_right, ROT_ANG);

%% 每条丝带对应的左右子弧（按占比累进切分）
segL = cell(1,K);
for k=1:K, segL{k} = subdivide_arc(theta_left(k,:),  F(:,k)); end    % 按列细分
segR = cell(1,P);
for i=1:P, segR{i} = subdivide_arc(theta_right(i,:), F(i,:)); end    % 按行细分

%% 作图
f = figure('Color','w','Position',[40 40 1100 820]); hold on;
ax = axes('Parent',f); hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');

% —— 左侧：簇（有厚度的节点弧带）+ 标签 ——
for k=1:K
    draw_arcband(theta_left(k,:), R_NODE_IN, R_NODE_OUT, Cclus(k,:), 0.95);
    place_label(theta_left(k,:), R_NODE_OUT + LABEL_OFFSET, sprintf('C%d',k), FONTSZ_LABEL, 'chord_label');
end

% —— 右侧：省份（有厚度的节点弧带）+ 标签（自动截断） ——
for i=1:P
    draw_arcband(theta_right(i,:), R_NODE_IN, R_NODE_OUT, Cprov(i,:), 0.95);
    txt_i = abbreviate_label(char(prov_top(i)), MAX_LABEL_LEN);
    place_label(theta_right(i,:), R_NODE_OUT + LABEL_OFFSET, txt_i, FONTSZ_LABEL, 'chord_label');
end

% —— 丝带（Province i → Cluster k；颜色按簇色） ——
for i=1:P
    for k=1:K
        if F(i,k)<=0, continue; end
        [thL, segL{k}] = pop_next(segL{k});
        [thR, segR{i}] = pop_next(segR{i});
        draw_ribbon(thL, thR, R_INNER, R_CTRL, Cclus(k,:), ALPHA_RIBBON);
    end
end

% —— 自动扩坐标范围以包含所有标签（防止溢出）——
lbl = findall(ax,'Type','text','Tag','chord_label');
if ~isempty(lbl)
    XY = [];
    for h = lbl.'
        set(h,'Units','data');
        ex = get(h,'Extent');                 % [x y w h] in data
        XY = [XY; ex(1), ex(2); ex(1)+ex(3), ex(2)+ex(4)]; %#ok<AGROW>
    end
    xmin = min(XY(:,1)); xmax = max(XY(:,1));
    ymin = min(XY(:,2)); ymax = max(XY(:,2));
    pad  = 0.05 * SCALE;   % 轻微留白
    xlim(ax,[xmin-pad, xmax+pad]); ylim(ax,[ymin-pad, ymax+pad]);
end

% —— 开启可视化编辑（可在图窗直接点选对象进行位置/文字编辑）——
plotedit on;
set(findall(ax,'Type','patch'),'PickableParts','all','HitTest','on');
set(findall(ax,'Type','text'), 'PickableParts','visible','HitTest','on');

% —— 可编辑标题（annotation，可鼠标双击编辑/拖动）——
if USE_WEIGHT, mode_str='weighted (sum w)'; else, mode_str='counts'; end
ht = annotation(f,'textbox',[0.15 0.93 0.70 0.06], ...
    'String', sprintf('Province → Cluster chord (Top-%d provinces, %s)', ...
                      min(TOPN_PROV,numel(prov_all)), mode_str), ...
    'EdgeColor','none','HorizontalAlignment','center', ...
    'FontWeight','bold','FontSize',12,'Interpreter','none','FitBoxToText','off');

% —— 保存 —— 

fname = fullfile(outdir, sprintf('fig2c3_chord_top%d_[%s]_%s', ...
         TOPN_PROV, cc_str, ternary_str(USE_WEIGHT,'weighted','counts')));

print(f,[fname,'.png'],'-dpng','-r300'); savefig(f,[fname,'.fig']);
disp(['Saved: ', fname, '.png & .fig']);

%% ==================== 辅助函数 ====================
function [labels, K] = get_labels_from_results(S)
    flds = fieldnames(S);
    if numel(flds)==1 && isstruct(S.(flds{1})), S=S.(flds{1}); end
    assert(isfield(S,'cluster') && isfield(S.cluster,'labels'), '缺少 cluster.labels');
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

function Tab = pick_widest_table(T)
    Tab=[]; fn=fieldnames(T); bestW=-1; key='';
    for i=1:numel(fn)
        if istable(T.(fn{i})) && width(T.(fn{i}))>bestW
            bestW=width(T.(fn{i})); key=fn{i};
        end
    end
    if ~isempty(key), Tab=T.(key); end
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

function s = tostr(x)
    if isstring(x), s=x; elseif ischar(x), s=string(x);
    elseif iscellstr(x), s=string(x);
    elseif iscell(x), try s=string(x); catch, s=strings(numel(x),1); for i=1:numel(x), if ischar(x{i}), s(i)=string(x{i}); end, end, end
    else, s=string(x); end
end

% —— 根据节点总量分配弧段（带组内间隙）
function th = alloc_arcs(weights, start_angle, span, gap)
    w = weights(:)'; w = max(w,0); S = sum(w);
    n = numel(w); th = zeros(n,2);
    if n==0, return; end
    if S==0
        base = span/max(n,1);
        for i=1:n
            a0 = start_angle + (i-1)*base;
            th(i,:) = [a0, a0 + base - gap];
        end
        return;
    end
    free_span = span - gap*(n-1);
    cur = start_angle;
    for i=1:n
        frac = w(i)/S;
        len  = frac * max(free_span,0);
        th(i,:) = [cur, cur+len];
        cur = cur + len + gap;
    end
end

% —— 按元素占比分割一段弧（返回队列 cell）
function Q = subdivide_arc(theta, parts)
    parts = parts(:)'; total = sum(parts);
    if total<=0, Q = {}; return; end
    L = theta(2)-theta(1);
    cur = theta(1);
    Q = cell(1,numel(parts));
    for j=1:numel(parts)
        len = (parts(j)/total)*L;
        Q{j} = [cur, cur+len];
        cur = cur + len;
    end
end

% —— 出队一个子弧
function [th, Qcell] = pop_next(Qcell)
    if isempty(Qcell) || isempty(Qcell{1}), th = []; return; end
    th = Qcell{1}; Qcell(1) = [];
end

% —— 画弧带（节点）
function draw_arcband(theta, r_in, r_out, col, faceAlpha)
    if nargin<5, faceAlpha=1; end
    t = linspace(theta(1), theta(2), 80);
    x1 = r_out*cos(t); y1 = r_out*sin(t);
    x2 = r_in *cos(fliplr(t)); y2 = r_in *sin(fliplr(t));
    patch([x1 x2], [y1 y2], col, 'EdgeColor','none', 'FaceAlpha', faceAlpha);
end

% —— 放标签（自动左右对齐） 
function place_label(theta, Rlab, txt, fsz, tagstr)
    if nargin<5, tagstr=''; end
    ang = mean(theta); x = Rlab*cos(ang); y = Rlab*sin(ang);
    if cos(ang) < 0, ha='right'; rot = rad2deg(ang)+180; else, ha='left'; rot = rad2deg(ang); end
    text(x, y, txt, 'Rotation', rot, 'HorizontalAlignment', ha, ...
         'VerticalAlignment','middle','FontSize',fsz,'Tag',tagstr);
end

% —— 画丝带（两端弧 thL / thR）
function draw_ribbon(thL, thR, r_in, r_ctrl, col, alpha)
    if isempty(thL) || isempty(thR), return; end
    tL0 = thL(1); tL1 = thL(2); tR0 = thR(1); tR1 = thR(2);
    [x1,y1] = pol2cart(tL0, r_in);  [x2,y2] = pol2cart(tR1, r_in);
    [x3,y3] = pol2cart(tL1, r_in);  [x4,y4] = pol2cart(tR0, r_in);
    [cx1,cy1] = pol2cart(tL0, r_ctrl); [cx2,cy2] = pol2cart(tR1, r_ctrl);
    [cx3,cy3] = pol2cart(tL1, r_ctrl); [cx4,cy4] = pol2cart(tR0, r_ctrl);
    N=30; tt=linspace(0,1,N);
    up   = bezier([x1 y1],[cx1 cy1],[cx2 cy2],[x2 y2], tt);
    down = bezier([x3 y3],[cx3 cy3],[cx4 cy4],[x4 y4], tt);
    poly = [up; flipud(down)];
    patch('XData',poly(:,1),'YData',poly(:,2), ...
          'FaceColor',col,'EdgeColor','none','FaceAlpha',alpha, ...
          'PickableParts','all','HitTest','on');
end

% —— 三次贝塞尔
function P = bezier(P0,P1,P2,P3,t)
    t = t(:); u=1-t;
    P = (u.^3).*P0 + 3*(u.^2).*t.*P1 + 3*u.*(t.^2).*P2 + (t.^3).*P3;
end

% —— 旋转角度矩阵
function TH2 = rot_theta(TH, ang)
    TH2 = TH + ang;
    TH2 = mod(TH2 + pi, 2*pi) - pi;  % 规范到 [-pi, pi)
end

% —— 省名截断
function txt2 = abbreviate_label(txt, maxlen)
    s = string(txt);
    if strlength(s) > maxlen, txt2 = char(extractBefore(s, maxlen+1) + "…");
    else, txt2 = char(s); end
end

% —— 三元字符串（无依赖）
function s = ternary_str(cond, a, b), if cond, s=a; else, s=b; end
end
