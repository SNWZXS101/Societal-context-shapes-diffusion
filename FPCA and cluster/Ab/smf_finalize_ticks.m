function handles = smf_finalize_ticks(target, varargin)
% smf_finalize_ticks  在编辑器里调完图后，一键“隐藏轴线并贴回刻度数字”
% 支持 figure 或 axes 作为 target；默认当前坐标轴 gca。
%
% 参数（Name-Value）：
%   'ShowX'        (true)     是否贴 X 刻度数字
%   'ShowY'        (true)     是否贴 Y 刻度数字
%   'Inside'       (false)    刻度数字放内侧(true)/外侧(false)
%   'OffsetFrac'   (0.02)     数字与边界的距离，占对应轴范围比例（0.01~0.05 常用）
%   'HideAxes'     (true)     是否把坐标轴隐藏（axis off）
%   'Bind'         (true)     数字字体与轴字体联动（改轴 FontSize/FontName 会同步）
%   'Color'        ([0 0 0])  刻度数字颜色
%   'UseTickLabels'(true)     使用当前 TickLabel（否则用数值转字符串）
%   'Remove'       (false)    仅清除手工贴的数字与监听器，不做其它操作
%
% 返回：
%   handles.X, handles.Y   两组 text 句柄（可能为空）

if nargin<1 || isempty(target), target = gca; end

% 解析参数
ip = inputParser;
ip.addParameter('ShowX', true);
ip.addParameter('ShowY', true);
ip.addParameter('Inside', false);
ip.addParameter('OffsetFrac', 0.02);
ip.addParameter('HideAxes', true);
ip.addParameter('Bind', true);
ip.addParameter('Color', [0 0 0]);
ip.addParameter('UseTickLabels', true);
ip.addParameter('Remove', false);
ip.parse(varargin{:});
op = ip.Results;

% 目标坐标轴集合
if isgraphics(target,'axes')
    axs = target;
elseif isgraphics(target,'figure')
    axs = findall(target,'Type','axes');  % 若需排除 colorbar/legend 可加条件
else
    error('target 必须是 axes 或 figure 句柄');
end
handles = struct('X',gobjects(0,1),'Y',gobjects(0,1));

for ax = reshape(axs,1,[])
    % 先清理旧的手工数字与监听器
    local_clear_manual(ax);

    if op.Remove
        continue; % 仅清除模式
    end

    % 需要时隐藏坐标轴（不再可点击）
    if op.HideAxes
        axis(ax, 'off');
    end

    % 读取范围/刻度
    xl = ax.XLim; yl = ax.YLim;
    xt = ax.XTick; yt = ax.YTick;

    % 刻度标签
    if op.UseTickLabels
        xtl = toStringVec(ax.XTickLabel, xt);
        ytl = toStringVec(ax.YTickLabel, yt);
    else
        xtl = string(xt);
        ytl = string(yt);
    end

    % 偏移与对齐
    dx = op.OffsetFrac * max(eps, (xl(2)-xl(1)));
    dy = op.OffsetFrac * max(eps, (yl(2)-yl(1)));
    if op.Inside
        yPosX = yl(1) + dy; vAlignX = 'bottom';
        xPosY = xl(1) + dx; hAlignY = 'left';
    else
        yPosX = yl(1) - dy; vAlignX = 'top';
        xPosY = xl(1) - dx; hAlignY = 'right';
    end

    % 默认字体取轴的
    fsize = ax.FontSize;
    fname = ax.FontName;

    hold(ax,'on');
    hX = gobjects(0); hY = gobjects(0);

    % 贴 X 数字
    if op.ShowX && ~isempty(xt)
        for i=1:numel(xt)
            if i<=numel(xtl) && xtl(i)~=""
                hX(end+1,1) = text(ax, xt(i), yPosX, 0, xtl(i), ...
                    'HorizontalAlignment','center', 'VerticalAlignment',vAlignX, ...
                    'FontName', fname, 'FontSize', fsize, 'Color', op.Color, ...
                    'Clipping','off', 'Tag','smfTickNumX');
            end
        end
    end

    % 贴 Y 数字
    if op.ShowY && ~isempty(yt)
        for i=1:numel(yt)
            if i<=numel(ytl) && ytl(i)~=""
                hY(end+1,1) = text(ax, xPosY, yt(i), 0, ytl(i), ...
                    'HorizontalAlignment', hAlignY, 'VerticalAlignment','middle', ...
                    'FontName', fname, 'FontSize', fsize, 'Color', op.Color, ...
                    'Clipping','off', 'Tag','smfTickNumY');
            end
        end
    end

    % 字体联动：监听轴字体变化
    if op.Bind
        set([hX;hY], 'FontSize', ax.FontSize, 'FontName', ax.FontName);
        L(1) = addlistener(ax,'FontSize','PostSet', @(~,~) set([hX;hY],'FontSize',ax.FontSize)); %#ok<AGROW>
        L(2) = addlistener(ax,'FontName','PostSet', @(~,~) set([hX;hY],'FontName',ax.FontName)); %#ok<AGROW>
        setappdata(ax, 'smfTickFontListeners', L);
    end

    % 累积返回
    handles.X = [handles.X; hX];
    handles.Y = [handles.Y; hY];
end
end

% --------- helpers ---------
function local_clear_manual(ax)
% 删除手工数字与旧监听器
delete(findall(ax,'Type','text','Tag','smfTickNumX'));
delete(findall(ax,'Type','text','Tag','smfTickNumY'));
L = getappdata(ax,'smfTickFontListeners');
if ~isempty(L)
    try, delete(L); end %#ok<TRYNC>
    rmappdata(ax,'smfTickFontListeners');
end
end

function s = toStringVec(lbl, ticks)
% 把可能是 char 矩阵/元胞/字符串的刻度标签转成 string 向量
if isstring(lbl)
    s = lbl;
elseif iscell(lbl)
    s = string(lbl);
elseif ischar(lbl)
    % char 矩阵：逐行取
    s = string(cellstr(lbl));
else
    s = string(ticks);
end
end
