function H = fig4a_cloud_rain_ppc(varargin)
% FIG4A_CLOUD_RAIN_PPC (fixed RGBA -> RGB+Alpha)
% Fig.4a Cloud–Rain PPC：左=每观测 P10–P90“雨线”+ y_obs 点；右= yhat_p50 vs y_obs 散点+拟合线
%
% 用法：
%   H = fig4a_cloud_rain_ppc('results_root','results_hbm', 'target','tf2_peak', ...
%         'export_pdf','Fig4a_tf2_peak.pdf');

% ---------- 参数 ----------
p = inputParser;
addParameter(p,'results_root','results_hbm',@(s)ischar(s)||isstring(s));
addParameter(p,'target','tf2_peak',@(s)ischar(s)||isstring(s));
addParameter(p,'export_pdf','',@(s)ischar(s)||isstring(s));
addParameter(p,'dpi',300,@(x)isnumeric(x)&&x>=72);
addParameter(p,'point_alpha',0.55,@(x)isnumeric(x)&&x>0&&x<=1);
addParameter(p,'line_alpha',0.75,@(x)isnumeric(x)&&x>0&&x<=1);
addParameter(p,'rain_color',[0 0 0],@(v)isnumeric(v)&&numel(v)==3);
addParameter(p,'obs_color',[0 0 0],@(v)isnumeric(v)&&numel(v)==3);
parse(p,varargin{:}); A = p.Results;

% ---------- 读取数据 ----------
ppc_path = fullfile(A.results_root,'M0_main',char(A.target),'ppc.csv');
if ~isfile(ppc_path), error('未找到 %s', ppc_path); end
T = readtable(ppc_path);
need = {'yhat_p10','yhat_p50','yhat_p90','yhat_mean','y_obs'};
for k=1:numel(need)
    if ~ismember(need{k}, T.Properties.VariableNames)
        error('ppc.csv 缺少列：%s', need{k});
    end
end
y10  = double(T.yhat_p10(:));
y50  = double(T.yhat_p50(:));
y90  = double(T.yhat_p90(:));
yhat = double(T.yhat_mean(:));
yobs = double(T.y_obs(:));
ok = isfinite(y10)&isfinite(y50)&isfinite(y90)&isfinite(yhat)&isfinite(yobs);
y10=y10(ok); y50=y50(ok); y90=y90(ok); yhat=yhat(ok); yobs=yobs(ok);
N = numel(yobs); if N==0, error('ppc.csv 有效行数为 0。'); end

% 排序（按 yhat_p50 从小到大）
[~, ord] = sort(y50,'ascend');
y10=y10(ord); y50=y50(ord); y90=y90(ord); yhat=yhat(ord); yobs=yobs(ord);

% 覆盖率
below  = mean(yobs < y10);
within = mean(yobs >= y10 & yobs <= y90);
above  = mean(yobs > y90);

% 右侧拟合与指标（在变换后尺度）
[slope, intercept, R2, MAE] = local_fit_stats(y50, yobs);

% ---------- 画布（自由布局，可拖拽） ----------
H.fig = figure('Color','w','Units','pixels','Position',[60 60 1280 620], ...
    'Name',sprintf('Fig4a Cloud–Rain PPC | %s', char(A.target)));
set(H.fig,'DefaultTextInterpreter','latex');
set(H.fig,'DefaultAxesTickLabelInterpreter','tex');
set(H.fig,'DefaultLegendInterpreter','latex');

H.ax(1) = axes('Parent',H.fig,'Units','normalized','Position',[0.08 0.14 0.64 0.78], ...
    'Box','off','LineWidth',1,'TickDir','out','FontName','Helvetica','Tag','ax_main');
H.ax(2) = axes('Parent',H.fig,'Units','normalized','Position',[0.77 0.14 0.19 0.78], ...
    'Box','off','LineWidth',1,'TickDir','out','FontName','Helvetica','Tag','ax_side');

plotedit(H.fig,'on'); safe_pick(H.ax(1)); safe_pick(H.ax(2));

% ---------- (a) 主面板 ----------
ax = H.ax(1); cla(ax); hold(ax,'on');
x = 1:N;

% 颜色：用“与白色混合”的方式模拟线条淡色（因为 line 不支持 alpha）
rainRGB  = blend_to_white(A.rain_color, 1 - A.line_alpha);   % 越靠近 1 越淡
midRGB   = blend_to_white([0 0 0], 0.75);                    % y50 的淡灰点
bandRGB  = blend_to_white([0 0 0], 0.85);                    % 分箱分位线更淡

% 雨线（P10–P90）
for i = 1:N
    plot(ax, [x(i) x(i)], [y10(i) y90(i)], '-', 'Color', rainRGB, 'LineWidth', 0.9, 'HitTest','on');
end
% yhat_p50 中点
plot(ax, x, y50, '.', 'Color', midRGB, 'MarkerSize', 12, 'HitTest','on');

% 观测点（透明度用 MarkerFaceAlpha）
scatter(ax, x, yobs, 9, 'filled', 'MarkerFaceColor', A.obs_color, ...
        'MarkerFaceAlpha', A.point_alpha, 'MarkerEdgeColor','none', 'HitTest','on');

% 局部分位带（x 轴等分箱内的 y_obs p10/p90）
nbin = max(12, round(sqrt(N)));
edges = round(linspace(1, N+1, nbin+1));
for b = 1:nbin
    idx = x>=edges(b) & x<edges(b+1);
    if ~any(idx), continue; end
    yy = yobs(idx); p10 = prctile(yy,10); p90 = prctile(yy,90);
    xb  = mean([edges(b), edges(b+1)]) - 0.5;
    dx  = 0.45*(edges(b+1)-edges(b));
    plot(ax, [xb-dx xb+dx], [p10 p10], '-', 'Color', bandRGB, 'LineWidth',1.2, 'HitTest','off');
    plot(ax, [xb-dx xb+dx], [p90 p90], '-', 'Color', bandRGB, 'LineWidth',1.2, 'HitTest','off');
end

% 轴与标题
xlim(ax, [0.5, N+0.5]);
yl = [min([y10; yobs])-eps, max([y90; yobs])+eps];
if yl(1)==yl(2), yl = yl + [-1 1]*0.5; end
ylim(ax, yl);
xlabel(ax, 'Observation (sorted by $\hat{y}_{50}$)','Interpreter','latex');
ylabel(ax, 'Value (transformed scale)');
title(ax, sprintf('Cloud–Rain PPC  \\it(%s)', char(A.target)), 'FontWeight','bold', 'FontSize',13);

% 覆盖率角注（TeX）
txt = sprintf(['90\\%% PI cover: within=\\textbf{%.1f}\\%% (exp. 90\\%%),  ', ...
               'below=%.1f\\%%, above=%.1f\\%%  |  N=%d'], ...
               within*100, below*100, above*100, N);
annotation(H.fig,'textbox',[0.08 0.93 0.84 0.05], 'String', txt, ...
    'EdgeColor','none','HorizontalAlignment','left','VerticalAlignment','middle', ...
    'FontName','Helvetica','FontSize',11,'Interpreter','latex');

% ---------- 右侧窄面板 ----------
ax2 = H.ax(2); cla(ax2); hold(ax2,'on');
scatter(ax2, y50, yobs, 28, 'filled', 'MarkerFaceColor', [0 0.1843 0.6588], ...
        'MarkerFaceAlpha', 0.35, 'MarkerEdgeColor','none', 'HitTest','on');

mm = [min([y50; yobs]), max([y50; yobs])];
if mm(1)==mm(2), mm = mm + [-1 1]*0.5; end
plot(ax2, mm, mm, '-', 'Color', blend_to_white([0 0 0],0.75), 'LineWidth', 1.3);

xx = linspace(mm(1), mm(2), 100);
yy = intercept + slope*xx;
plot(ax2, xx, yy, '-', 'Color', [0 0 0], 'LineWidth', 1.4);

axis(ax2,'tight'); xlim(ax2, mm); ylim(ax2, mm);
xlabel(ax2, '$\hat{y}_{50}$','Interpreter','latex');
ylabel(ax2, 'y_{\rm obs}');
title(ax2, 'Fit: slope/intercept, R^2, MAE', 'FontWeight','bold','FontSize',13,'Interpreter','tex');

txt2 = sprintf('slope=%.3f,  intercept=%.3f\nR^2=%.3f,  MAE=%.3f', slope, intercept, R2, MAE);
annotation(H.fig,'textbox',[0.77 0.92 0.21 0.07],'String',txt2, ...
    'EdgeColor','none','HorizontalAlignment','left','VerticalAlignment','middle', ...
    'FontName','Helvetica','FontSize',11,'Interpreter','tex');

set([ax,ax2],'TickLabelInterpreter','tex','FontName','Helvetica');
grid(ax,'off'); grid(ax2,'off');

% ---------- 导出（可选） ----------
if ~isempty(A.export_pdf)
    try
        set(H.fig,'Renderer','painters','PaperPositionMode','auto');
        print(H.fig, A.export_pdf, '-dpdf', '-painters', '-bestfit');
        fprintf('[Fig4a] 导出 PDF：%s\n', A.export_pdf);
    catch ME
        warning(sprintf('导出失败：%s', ME.message));
    end
end
end

% ======= 工具函数 =======
function [slope, intercept, R2, MAE] = local_fit_stats(x, y)
X = [ones(numel(x),1), x(:)];
b = X\y(:);
intercept = b(1); slope = b(2);
yfit = X*b;
res  = y(:) - yfit;
SSres = sum(res.^2);
SStot = sum( (y(:) - mean(y(:))).^2 );
R2 = 1 - SSres/SStot;
MAE = mean(abs(res));
end

function c = blend_to_white(rgb, t)
% 把颜色 rgb 与白色按 t 混合（t∈[0,1]：0=原色；1=全白）
t = min(max(t,0),1);
c = (1-t)*rgb + t*[1 1 1];
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
