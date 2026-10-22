function S = make_panels_editable_for_timeline(H, varargin)
% MAKE_PANELS_EDITABLE 让 tiledlayout 里的各子图可自由拖动 / 改大小
%
% 用法：
%   H = struct('fig', gcf, 'ax', findall(gcf,'Type','axes'));
%   make_panels_editable(H);
%
%   % 或直接：
%   make_panels_editable(gcf);              % 传 figure 句柄
%   make_panels_editable(findall(gcf,'Type','axes'));  % 传 axes 数组
%
% 可选参数（Name-Value）：
%   'DeleteTiled'  (true)  移除找到的 tiledlayout 容器（推荐）
%   'KeepOrder'    (true)  维持当前视觉顺序（从左上到右下）
%   'PickableAll'  (true)  打开所有对象的 HitTest/PickableParts
%   'OpenEditor'   (true)  结束时打开 plotedit 编辑模式
%
% 返回：
%   S：结构体，保存了原对象与位置，便于需要时恢复（S.restore()）
%
% 注意：
%   - 该函数会把每个 axes 的 Parent 改为 figure、自适应 normalized 坐标；
%   - legend / colorbar 若是“轴附属”对象，会一并移动并保持相对位置；
%   - 页眉（sgtitle）是 annotation/text 容器，不会阻止拖动，必要时可手动删。
%
% 作者：你的 MATLAB 可视化小助手（2025）

% -------- 参数解析 --------
if nargin == 0 || isempty(H)
    H = struct('fig', gcf, 'ax', findall(gcf,'Type','axes'));
end

p = inputParser;
addParameter(p,'DeleteTiled', true,  @(x)islogical(x)&&isscalar(x));
addParameter(p,'KeepOrder',   true,  @(x)islogical(x)&&isscalar(x));
addParameter(p,'PickableAll', true,  @(x)islogical(x)&&isscalar(x));
addParameter(p,'OpenEditor',  true,  @(x)islogical(x)&&isscalar(x));
parse(p, varargin{:});
opt = p.Results;

% -------- 规格化输入 --------
if ishghandle(H)
    if strcmpi(get(H,'Type'),'figure')
        fig = H;
        ax  = findall(fig,'Type','axes');
    elseif strcmpi(get(H,'Type'),'axes')
        ax  = H(:);
        fig = ancestor(ax(1),'figure');
    else
        error('不支持的句柄类型：%s', get(H,'Type'));
    end
elseif isstruct(H) && isfield(H,'fig') && isfield(H,'ax')
    fig = H.fig;
    ax  = H.ax(:);
else
    error('输入应为 figure 句柄、axes 数组，或含 fig/ax 的结构体。');
end

if isempty(fig) || ~ishandle(fig), error('H.fig 无效'); end
ax = ax(ishandle(ax));
if isempty(ax), error('未找到有效的 axes。'); end

% 可选：按当前视觉顺序（从上到下、左到右）排序，复制时更直观
if opt.KeepOrder
    ax = sort_axes_by_visual_order(ax);
end

% 顶层 figure 像素尺寸
pf = getpixelposition(fig, true);
if pf(3)<=0 || pf(4)<=0
    drawnow; pf = getpixelposition(fig, true);
end

% -------- 收集所有相关对象（axes + legend/colorbar）--------
N = numel(ax);
AxInfo = struct('ax',[],'pos_px',[],'parent_before',[],'legend',[],'cbar',[],'tlo',[]);
AxInfo = repmat(AxInfo, N, 1);

% 找可能的父级 tiledlayout（每个轴各自可能不同）
for k = 1:N
    AxInfo(k).ax = ax(k);
    AxInfo(k).parent_before = ax(k).Parent;
    AxInfo(k).pos_px = getpixelposition(ax(k), true);

    % 关联 legend（若有）
    lg = legend(ax(k));  % 仅获取已存在的 legend 句柄
    if ~isempty(lg) && isvalid(lg) && strcmp(lg.Visible,'on')
        AxInfo(k).legend = lg;
        AxInfo(k).legend_pos_px = getpixelposition(lg, true);
    end

    % 关联 colorbar（若有）
    % 注意：可能存在多个 colorbar，这里只取“与该轴绑定”的
    cb = findall(fig,'Type','ColorBar','-and','Axes',ax(k));
    if ~isempty(cb)
        cb = cb(1);
        AxInfo(k).cbar = cb;
        AxInfo(k).cbar_pos_px = getpixelposition(cb, true);
    end

    % 追溯 tiledlayout
    try
        tlo = ancestor(ax(k), 'matlab.graphics.layout.TiledChartLayout');
        if ~isempty(tlo) && isvalid(tlo)
            AxInfo(k).tlo = tlo;
        end
    catch
        AxInfo(k).tlo = [];
    end
end

% -------- 位置与父子关系迁移（核心步骤）--------
for k = 1:N
    axk = AxInfo(k).ax;
    pk  = AxInfo(k).pos_px;

    % 1) 归一化坐标
    posN = [pk(1)/pf(3), pk(2)/pf(4), pk(3)/pf(3), pk(4)/pf(4)];
    set(axk, 'Units','normalized', 'Position',posN);

    % 2) 重新挂到 figure（脱离 tiledlayout）
    axk.Parent = fig;

    % 3) legend / colorbar 同步位置（若存在）
    if ~isempty(AxInfo(k).legend) && isvalid(AxInfo(k).legend)
        pl = AxInfo(k).legend_pos_px;
        posNl = [pl(1)/pf(3), pl(2)/pf(4), pl(3)/pf(3), pl(4)/pf(4)];
        set(AxInfo(k).legend, 'Units','normalized', 'Position',posNl);
        AxInfo(k).legend.Parent = fig;
    end
    if ~isempty(AxInfo(k).cbar) && isvalid(AxInfo(k).cbar)
        pc = AxInfo(k).cbar_pos_px;
        posNc = [pc(1)/pf(3), pc(2)/pf(4), pc(3)/pf(3), pc(4)/pf(4)];
        set(AxInfo(k).cbar, 'Units','normalized', 'Position',posNc);
        AxInfo(k).cbar.Parent = fig;
    end
end

% -------- 可选：删除 tiledlayout 容器 --------
if opt.DeleteTiled
    % 收集唯一的布局容器并删除
    TL = unique([AxInfo.tlo]);
    TL = TL(isgraphics(TL));   % 用基础图形有效性判断，避免 isvalid 冲突
    for i = 1:numel(TL)
        try
            delete(TL(i));
        catch
            % 忽略删除失败的情况
        end
    end
    for i = 1:numel(TL)
        try, delete(TL(i)); catch, end
    end
end

% -------- 打开选择/编辑能力 --------
set(fig, 'Renderer','painters');  % 确保矢量可编辑
if opt.PickableAll
    kids = findall(fig,'-property','HitTest');
    set(kids, 'HitTest','on');
    kids = findall(fig,'-property','PickableParts');
    set(kids, 'PickableParts','all');
end
if opt.OpenEditor
    try, plotedit(fig,'on'); catch, end
end

% -------- 返回可恢复信息 --------
S = struct();
S.fig = fig;
S.axes = ax;
S.info = AxInfo;
S.restore = @()restore_layout(S);

end % ===== 主函数结束 =====


% ----------------- 辅助：恢复到原父对象（尽力而为） -----------------
function restore_layout(S)
if ~isstruct(S) || ~isfield(S,'info'), return; end
try
    fig = S.fig;
    pf  = getpixelposition(fig, true);
    for k = 1:numel(S.info)
        if ~ishandle(S.info(k).ax), continue; end
        axk = S.info(k).ax;
        % 恢复到原父级（若仍存在）
        if isvalid(S.info(k).parent_before)
            axk.Parent = S.info(k).parent_before;
        else
            axk.Parent = fig;
        end
        % 尽力恢复像素位置（父容器不同，可能略有偏差）
        pk = S.info(k).pos_px;
        if all(pk(3:4) > 0)
            setpixelposition(axk, pk, true);
        end
        % legend / colorbar
        if ~isempty(S.info(k).legend) && isvalid(S.info(k).legend)
            setpixelposition(S.info(k).legend, S.info(k).legend_pos_px, true);
            if isvalid(S.info(k).parent_before)
                S.info(k).legend.Parent = S.info(k).parent_before;
            else
                S.info(k).legend.Parent = fig;
            end
        end
        if ~isempty(S.info(k).cbar) && isvalid(S.info(k).cbar)
            setpixelposition(S.info(k).cbar, S.info(k).cbar_pos_px, true);
            if isvalid(S.info(k).parent_before)
                S.info(k).cbar.Parent = S.info(k).parent_before;
            else
                S.info(k).cbar.Parent = fig;
            end
        end
    end
catch
    % 忽略恢复中的非致命问题
end
end


% ----------------- 辅助：按视觉顺序排序 axes -----------------
function ax_sorted = sort_axes_by_visual_order(ax)
ax = ax(:); n = numel(ax);
pos = zeros(n,4);
for i=1:n
    try
        pos(i,:) = getpixelposition(ax(i), true);
    catch
        pos(i,:) = [0 0 0 0];
    end
end
% 视觉上：先按 y 从大到小（顶到下），再按 x 从小到大（左到右）
[~, idx] = sortrows([-(pos(:,2)+pos(:,4)) , pos(:,1)]);
ax_sorted = ax(idx);
end
