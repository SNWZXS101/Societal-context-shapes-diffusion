function H = fig4b_residual_raincloud_by_variant(varargin)
% FIG4B_RESIDUAL_RAINCLOUD_BY_VARIANT (bw fixed)
% 残差雨云图（半小提琴 + 箱线 + 抖动），按变异株分组：
%   r = y_obs - yhat_p50（或 yhat_mean），全在“变换后尺度”。

% ---------- 参数 ----------
p = inputParser;
addParameter(p,'results_root','results_hbm',@(s)ischar(s)||isstring(s));
addParameter(p,'target','tf2_peak',@(s)ischar(s)||isstring(s));
addParameter(p,'resid_from','median',@(s)any(strcmpi(s,{'median','mean'})));
addParameter(p,'export_pdf','',@(s)ischar(s)||isstring(s));
addParameter(p,'violin_width',0.34,@(x)isnumeric(x)&&x>0);
addParameter(p,'point_alpha',0.60,@(x)isnumeric(x)&&x>0&&x<=1);
addParameter(p,'point_size',25,@(x)isnumeric(x)&&x>0);
addParameter(p,'jitter',0.06,@(x)isnumeric(x)&&x>=0);
addParameter(p,'bw','plug-in',@(x)ischar(x)||isstring(x)||isnumeric(x)); % ← 默认 plug-in
parse(p,varargin{:}); A = p.Results;

% ---------- 读数据 ----------
ppc_path   = fullfile(A.results_root,'M0_main',char(A.target),'ppc.csv');
cache_path = fullfile(A.results_root,'cache','hbm_cache.mat');
if ~isfile(ppc_path),   error('未找到 %s', ppc_path); end
if ~isfile(cache_path), error('未找到 %s', cache_path); end
Tppc  = readtable(ppc_path);
S     = load(cache_path); cache = S.cache;

y_tr  = cache.Y.(char(A.target));
ok    = isfinite(y_tr) & all(isfinite(cache.X_fe),2);
var_use = cache.var_id(ok);
var_names = string(cache.var_levels(:));
N_use = sum(ok);
if height(Tppc) ~= N_use
    error('ppc.csv 行数(%d)与可用行数(%d)不一致；请确认导出版本。', height(Tppc), N_use);
end

switch lower(A.resid_from)
    case 'median', yhat_ref = Tppc.yhat_p50;
    otherwise,     yhat_ref = Tppc.yhat_mean;
end
yobs = Tppc.y_obs;
r    = double(yobs - yhat_ref);

want = intersect(["Original","Delta","Omicron"], var_names, 'stable');
if isempty(want), want = var_names(:)'; end
G = want;

data = cell(1,numel(G));
for i=1:numel(G)
    gid  = find(var_names==G(i));
    data{i} = r(var_use==gid);
end

% ---------- 颜色 ----------
col_original = [86 140 199]/255;   % 蓝
col_delta    = [233 135  59]/255;  % 橙
col_omicron  = [111 191 115]/255;  % 绿
Cmap = containers.Map({'Original','Delta','Omicron'}, ...
                      {col_original,col_delta,col_omicron});
cols = zeros(numel(G),3);
for i=1:numel(G)
    k = char(G(i));
    if isKey(Cmap,k), cols(i,:)=Cmap(k); else, cols(i,:)=[0.4 0.4 0.4]; end
end

% ---------- 画布（可拖拽） ----------
H.fig = figure('Color','w','Units','pixels','Position',[80 80 1020 520], ...
    'Name',sprintf('Fig4b Residual raincloud | %s', char(A.target)));
set(H.fig,'DefaultTextInterpreter','tex');
set(H.fig,'DefaultAxesTickLabelInterpreter','tex');
set(H.fig,'DefaultLegendInterpreter','tex');

H.ax(1) = axes('Parent',H.fig,'Units','normalized','Position',[0.09 0.16 0.86 0.78], ...
    'Box','off','LineWidth',1.5,'TickDir','out','FontName','Helvetica','Tag','ax_main','FontSize',15);
plotedit(H.fig,'on'); safe_pick(H.ax(1));

ax = H.ax(1); cla(ax); hold(ax,'on');

% 0 线
plot(ax, xlim(ax), [0 0], '-', 'Color',[0 0 0 0.18], 'LineWidth',1.2, 'HitTest','off');

% ---------- 带宽选项映射 ----------
[bwFlag, bwArg] = map_bandwidth(A.bw);

% ---------- 雨云 ----------
xpos = 1:numel(G);
for i=1:numel(G)
    xi = xpos(i);
    ri = data{i};
    ri = ri(isfinite(ri));
    if isempty(ri), continue; end

    % 半小提琴（左半）
    if numel(ri) >= 10
        if bwFlag==0      % 默认/自动
            [yy, ff] = ksdensity(ri);
        else              % 指定带宽或方法
            [yy, ff] = ksdensity(ri, bwArg{:});
        end
    else
        % 小样本兜底：正态核，带宽按 0.9*min(sd, IQR/1.34)*n^(-1/5)
        s = std(ri);
        iqr = iqr_robust(ri);
        h = 0.9*min(s, iqr/1.34) * max(numel(ri),1)^(-1/5);
        h = max(h, 1e-6);
        [yy, ff] = ksdensity(ri, 'Bandwidth', h);
    end
    ff = ff / max(ff + eps);
    vv = A.violin_width * ff;

    patch('Parent',ax, ...
          'XData', xi - vv, 'YData', yy, ...
          'FaceColor', cols(i,:), 'FaceAlpha', 0.62, ...
          'EdgeColor', 'none', 'HitTest','on');

    % 箱线
    q  = quantile(ri,[0.25 0.5 0.75]); iqrV = q(3)-q(1);
    lo = max(min(ri), q(1) - 1.5*iqrV);
    hi = min(max(ri), q(3) + 1.5*iqrV);
    rectangle('Position',[xi-0.10 q(1) 0.20 (q(3)-q(1))], ...
              'FaceColor',[1 1 1], 'EdgeColor',[0 0 0], 'LineWidth',0.8, 'Parent',ax);
    plot(ax, [xi-0.10 xi+0.10], [q(2) q(2)], '-', 'Color',[0 0 0], 'LineWidth',1.2);
    plot(ax, [xi xi], [q(3) hi], '-', 'Color',[0 0 0], 'LineWidth',0.8);
    plot(ax, [xi xi], [q(1) lo], '-', 'Color',[0 0 0], 'LineWidth',0.8);
    plot(ax, [xi-0.05 xi+0.05], [hi hi], '-', 'Color',[0 0 0], 'LineWidth',0.8);
    plot(ax, [xi-0.05 xi+0.05], [lo lo], '-', 'Color',[0 0 0], 'LineWidth',0.8);

    % 抖动点（右侧）
    n = numel(ri);
    jx = xi + 0.12 + A.jitter*(2*rand(n,1)-1);
    scatter(ax, jx, ri, A.point_size, 'filled', ...
            'MarkerFaceColor', cols(i,:), 'MarkerFaceAlpha', A.point_alpha, ...
            'MarkerEdgeColor','none', 'HitTest','on');
end

% 轴 & 文本
xlim(ax, [0.5, numel(G)+0.5]);
allr = cell2mat(data(:)); if isempty(allr), allr = 0; end
yr = prctile(allr, [1 99]); d = diff(yr);
ylim(ax, [yr(1)-0.10*d, yr(2)+0.10*d]);

set(ax,'XTick',xpos,'XTickLabel',strrep(cellstr(G),'_','\_'));
set(ax,'YGrid','on'); ax.GridAlpha = 0.08;

% y 标签（含 \hat，单独 latex）
if strcmpi(A.resid_from,'median')
    ylabel(ax, '$r = y_{\mathrm{obs}} - \hat{y}_{50}$', 'Interpreter','latex');
else
    ylabel(ax, '$r = y_{\mathrm{obs}} - \bar{y}_{\mathrm{rep}}$', 'Interpreter','latex');
end
xlabel(ax, 'Variant', 'Interpreter','tex');

title(ax, sprintf('Residual raincloud  \\it(%s)', char(A.target)), ...
      'FontWeight','bold','FontSize',16,'Interpreter','tex');

% 统计提示
grp = arrayfun(@(i) i*ones(numel(data{i}),1), 1:numel(G), 'uni',0);
grp = cat(1, grp{:}); vals = cat(1, data{:});
try
    p = kruskalwallis(vals, grp, 'off');
    txt = sprintf('Kruskal–Wallis:  \\itp\\rm=%.3g', p);
catch
    txt = 'Kruskal–Wallis:  n/a';
end
annotation(H.fig,'textbox',[0.09 0.93 0.82 0.05],'String',txt, ...
    'EdgeColor','none','HorizontalAlignment','left','VerticalAlignment','middle', ...
    'FontName','Helvetica','FontSize',12,'Interpreter','tex');

set(ax,'TickLabelInterpreter','tex','FontName','Helvetica','LineWidth',1.2,'TickDir','out');
box(ax,'off');

% 导出
if ~isempty(A.export_pdf)
    try
        set(H.fig,'Renderer','painters','PaperPositionMode','auto');
        print(H.fig, A.export_pdf, '-dpdf', '-painters', '-bestfit');
        fprintf('[Fig4b] 导出 PDF：%s\n', A.export_pdf);
    catch ME
        warning(sprintf('导出失败：%s', ME.message));
    end
end
end

% ===== 工具 =====
function [flag, arg] = map_bandwidth(bw)
% flag: 0=默认, 1=指定
if isnumeric(bw)
    flag = 1; arg = {'Bandwidth', max(bw,1e-9)};
    return;
end
bw = lower(string(bw));
switch bw
    case "auto"           % 兼容旧调用：映射为 plug-in
        flag = 1; arg = {'Bandwidth','plug-in'};
    case {"plug-in","plugin","plug in"}
        flag = 1; arg = {'Bandwidth','plug-in'};
    case {"normal-approx","normal","silverman"}
        flag = 1; arg = {'Bandwidth','normal-approx'};
    otherwise
        % 未知字符串：回退到默认（由 ksdensity 自定）
        flag = 0; arg = {};
end
end

function v = iqr_robust(x)
q = quantile(x,[0.25 0.75]); v = q(2)-q(1);
if ~isfinite(v) || v<=0, v = std(x); end
if ~isfinite(v) || v<=0, v = 1; end
end

function safe_pick(h)
try
    if isgraphics(h)
        if isprop(h,'HitTest'), set(h,'HitTest','on'); end
        if isprop(h,'PickableParts'), set(h,'PickableParts','all'); end
    end
catch
end
end
