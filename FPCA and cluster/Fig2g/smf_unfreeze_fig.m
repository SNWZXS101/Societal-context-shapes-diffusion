function smf_unfreeze_fig(hf)
% 解除 R2022b 交互卡顿：关闭交互、让 patch 不抢点击、轴可选中
if nargin<1 || isempty(hf), hf = gcf; end
try, datacursormode(hf,'off'); catch, end
try, brush(hf,'off'); catch, end
try, rotate3d(hf,'off'); catch, end
try, plotedit(hf,'off'); catch, end
try, zoom(hf,'off'); pan(hf,'off'); catch, end
set(hf,'WindowButtonDownFcn',[]);
set(hf,'WindowButtonUpFcn',[]);
set(hf,'WindowButtonMotionFcn',[]);
axs   = findall(hf,'Type','axes');
fills = findall(hf,'Type','patch');
set(fills,'HitTest','off','PickableParts','none');
set(axs,'HitTest','on','PickableParts','all');
try, arrayfun(@(ax) disableDefaultInteractivity(ax), axs); catch, end
try, set(axs,'Toolbar',[]); catch, end
disp('✅ smf_unfreeze_fig: 已解除交互卡顿。');
end
