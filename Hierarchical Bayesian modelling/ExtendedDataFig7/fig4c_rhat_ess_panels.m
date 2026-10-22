function H = fig4c_rhat_ess_panels(varargin)
% FIG4C_RHAT_ESS_PANELS
% Fig.4c: 采样诊断（R̂ & ESS）双面板，可拖拽编辑、TeX 文本。
%
% 用法示例：
%   H = fig4c_rhat_ess_panels('results_root','results_hbm', ...
%         'target','tf2_peak', 'export_pdf','Fig4c_tf2_peak.pdf');

% ---------- 参数 ----------
p = inputParser;
addParameter(p,'results_root','results_hbm',@(s)ischar(s)||isstring(s));
addParameter(p,'target','tf2_peak',@(s)ischar(s)||isstring(s));
addParameter(p,'export_pdf','',@(s)ischar(s)||isstring(s));
addParameter(p,'rhat_xlim',[0.995 1.08],@(v)isnumeric(v)&&numel(v)==2);
addParameter(p,'ess_log',false,@islogical);          % 如 ESS 范围极宽可开对数
addParameter(p,'ess_bins',30,@(x)isnumeric(x)&&x>=5);
addParameter(p,'rhat_bins',30,@(x)isnumeric(x)&&x>=5);
parse(p,varargin{:}); A = p.Results;

% ---------- 读数据 ----------
sum_path = fullfile(A.results_root,'M0_main',char(A.target),'summary.csv');
if ~isfile(sum_path)
    error('未找到 %s', sum_path);
end
T = readtable(sum_path);

% 容错列名
col_rhat = pick_col(T, {'rhat','Rhat','R_hat','split_rhat'});
col_ess  = pick_col(T, {'ess','n_eff','ess_bulk','ess_tail'});

if isempty(col_rhat) || isempty(col_ess)
    error('summary.csv 未找到 rhat/ess 列（已尝试多种候选名）。');
end

rhat = double(T.(col_rhat));
ess  = double(T.(col_ess));

% 清洗
rhat = rhat(isfinite(rhat) & rhat>0);
ess  = ess(isfinite(ess)   & ess>0);

if isempty(rhat) || isempty(ess)
    error('rhat 或 ess 有效数据为空，请检查 summary.csv。');
end

% ---------- 统计量 ----------
pct_lt_101 = mean(rhat < 1.01)*100;
pct_lt_105 = mean(rhat < 1.05)*100;
rhat_med   = median(rhat,'omitnan');
rhat_max   = max(rhat);

ess_med = median(ess,'omitnan');
ess_min = min(ess);
ess_q10 = quantile(ess,0.10);
ess_q25 = quantile(ess,0.25);

% ---------- 画布（自由布局，可拖拽） ----------
H.fig = figure('Color','w','Units','pixels','Position',[80 80 1180 520], ...
    'Name',sprintf('Fig4c R̂ & ESS | %s', char(A.target)));

% 全局 TeX
set(H.fig,'DefaultTextInterpreter','tex');
set(H.fig,'DefaultAxesTickLabelInterpreter','tex');
set(H.fig,'DefaultLegendInterpreter','tex');

% 两块 axes：左 Rhat，右 ESS
H.ax(1) = axes('Parent',H.fig,'Units','normalized','Position',[0.08 0.15 0.40 0.78], ...
    'Box','off','LineWidth',0.8,'TickDir','out','FontName','Helvetica','Tag','ax_rhat','FontSize',11);
H.ax(2) = axes('Parent',H.fig,'Units','normalized','Position',[0.56 0.15 0.40 0.78], ...
    'Box','off','LineWidth',0.8,'TickDir','out','FontName','Helvetica','Tag','ax_ess','FontSize',11);

plotedit(H.fig,'on'); safe_pick(H.ax(1)); safe_pick(H.ax(2));

% ---------- 左：R̂ 直方图 ----------
ax1 = H.ax(1); cla(ax1); hold(ax1,'on');

% 夹取到可视范围，避免极端值拉伸
xlim1 = A.rhat_xlim;
rhat_plot = rhat(rhat>=xlim1(1) & rhat<=xlim1(2));

% 直方图（柔灰）
h1 = histogram(ax1, rhat_plot, A.rhat_bins, ...
    'FaceColor',[0.35 0.45 0.60], 'FaceAlpha',0.55, 'EdgeColor','none', 'DisplayStyle','bar');
set(h1,'HitTest','on','PickableParts','all');

% 阈值线
plot(ax1, [1.01 1.01], ylim(ax1), '-', 'Color',[0 0 0 0.35], 'LineWidth',1.0, 'HitTest','off');
plot(ax1, [1.05 1.05], ylim(ax1), '--', 'Color',[0 0 0 0.35], 'LineWidth',1.0, 'HitTest','off');

xlim(ax1, xlim1);
xlabel(ax1, 'R̂','Interpreter','tex');   % R + U+0302
ylabel(ax1, 'Count','Interpreter','tex');
title(ax1, sprintf('R̂: <1.01 = %.1f%%,  <1.05 = %.1f%%', pct_lt_101, pct_lt_105), ...
      'FontWeight','bold','FontSize',12,'Interpreter','tex');


% 角注：中位/最大
txt1 = sprintf('median=%.3f   max=%.3f   n=%d', rhat_med, rhat_max, numel(rhat));
annotation(H.fig,'textbox',[0.08 0.93 0.40 0.06],'String',txt1, ...
    'EdgeColor','none','HorizontalAlignment','left','VerticalAlignment','middle', ...
    'FontName','Helvetica','FontSize',11,'Interpreter','tex');

% ---------- 右：ESS 直方图 ----------
ax2 = H.ax(2); cla(ax2); hold(ax2,'on');

ess_plot = ess;
if A.ess_log
    % 对数刻度：横轴对数，bin 采用对数等距边界
    set(ax2,'XScale','log');
    % 自动边界（避免 0）
    xmin = max(min(ess_plot), eps);
    xmax = max(ess_plot);
    edges = logspace(log10(xmin), log10(xmax), A.ess_bins+1);
    h2 = histogram(ax2, ess_plot, edges, ...
        'FaceColor',[0.30 0.55 0.45], 'FaceAlpha',0.55, 'EdgeColor','none');
else
    h2 = histogram(ax2, ess_plot, A.ess_bins, ...
        'FaceColor',[0.30 0.55 0.45], 'FaceAlpha',0.55, 'EdgeColor','none');
end
set(h2,'HitTest','on','PickableParts','all');

xlabel(ax2, 'ESS','Interpreter','tex');
ylabel(ax2, 'Count','Interpreter','tex');
title(ax2, sprintf('ESS: median=%.0f,  min=%.0f  (Q_{10}=%.0f, Q_{25}=%.0f)', ...
      ess_med, ess_min, ess_q10, ess_q25), ...
      'FontWeight','bold','FontSize',12,'Interpreter','tex');

% 微网格 & 去上右框
set([ax1,ax2],'YGrid','on'); ax1.GridAlpha=0.08; ax2.GridAlpha=0.08;
box(ax1,'off'); box(ax2,'off');

% 顶部总标题（可选）
sg = sgtitle(H.fig, sprintf('Sampler diagnostics | %s', char(A.target)));
set(sg,'FontSize',14,'FontWeight','bold','Interpreter','tex'); safe_pick(sg);

% ---------- 导出（可选） ----------
if ~isempty(A.export_pdf)
    try
        set(H.fig,'Renderer','painters','PaperPositionMode','auto');
        print(H.fig, A.export_pdf, '-dpdf', '-painters', '-bestfit');
        fprintf('[Fig4c] 导出 PDF：%s\n', A.export_pdf);
    catch ME
        warning(springf('导出失败：%s', ME.message));
    end
end
end

% ===== 工具 =====
function name = pick_col(T, cands)
name = '';
vn = string(T.Properties.VariableNames);
for i=1:numel(cands)
    hit = find(strcmpi(vn, cands{i}), 1);
    if ~isempty(hit), name = char(vn(hit)); return; end
end
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
