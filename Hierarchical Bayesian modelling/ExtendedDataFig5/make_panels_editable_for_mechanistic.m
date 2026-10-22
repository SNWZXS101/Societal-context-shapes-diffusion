function make_panels_editable_for_mechanistic(H)
% 让 tiledlayout 里的各子图变成可自由拖动/改大小的轴
% 用法：
%   H = plot_mechanistic_ridge_minmax1(...);
%   make_panels_editable(H);

fig = H.fig;
ax  = H.ax;
if isempty(fig) || ~ishandle(fig), error('H.fig 无效'); end
if isempty(ax)  || ~all(ishandle(ax(:))), error('H.ax 无效'); end

% 1) 找到顶层 tiledlayout（如果有）
tlo = [];
try
    tlo = ancestor(ax(1), 'matlab.graphics.layout.TiledChartLayout');
catch
    % 旧版/已脱离时可能找不到
end

% 2) 把每个轴的位置从像素坐标换算为 figure 归一化坐标，并重新挂到 figure
pf = getpixelposition(fig,true);
for k = 1:numel(ax)
    if ~ishandle(ax(k)), continue; end
    % 取当前“视觉上的”像素位置
    pk = getpixelposition(ax(k), true);
    % 归一化到 figure
    posN = [pk(1)/pf(3), pk(2)/pf(4), pk(3)/pf(3), pk(4)/pf(4)];
    ax(k).Units    = 'normalized';
    ax(k).Position = posN;
    ax(k).Parent   = fig;      % —— 脱离 tiledlayout 的关键
end

% 3) 允许交互选择/编辑（面对象、线对象、文字都可点选）
set(fig, 'Renderer', 'painters');     % 确保矢量可编辑
set(findall(fig,'-property','PickableParts'), 'PickableParts','all');
set(findall(fig,'-property','HitTest'),       'HitTest','on');

% 4) 删除（或隐藏）tiledlayout 容器，避免它继续管理布局
if ~isempty(tlo) && isvalid(tlo)
    delete(tlo);
end

% 5) 开启图形编辑模式（可用鼠标拖动/缩放对象和坐标轴）
try
    plotedit(fig,'on');   % 老指令，编辑器开关
catch
    % R2024b 也支持 plotedit；若遇到 UI 问题可考虑：plottools(fig)
end

end
