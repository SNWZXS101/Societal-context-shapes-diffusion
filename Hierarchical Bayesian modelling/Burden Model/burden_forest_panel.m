function burden_forest_panel(results_root, target_dir, panel_title, varargin)
% BURDEN_FOREST_PANEL
% 固定效应（beta）森林图：相对变化(%) = (exp(beta) - 1) * 100

p=inputParser;
addParameter(p,'Kshow',25,@(x)isnumeric(x)&&x>0);
addParameter(p,'Color',[0.12 0.45 0.75],@(x)isnumeric(x)&&numel(x)==3);
addParameter(p,'Out','',@(s)ischar(s)||isstring(s));
addParameter(p,'FontName','Arial',@(s)ischar(s)||isstring(s));
addParameter(p,'DPI',400,@(x)isnumeric(x)&&x>0);
parse(p,varargin{:});
opt=p.Results;

tbl = readtable(fullfile(results_root,'burden',target_dir,'summary.csv'));
ix  = startsWith(string(tbl.param),'beta:');
if ~any(ix), error('summary.csv 未找到 beta:* 行'); end
T   = tbl(ix,:);

% 相对变化（%）
m   = 100*(exp(T.mean)-1);
lo  = 100*(exp(T.ci2_5)-1);
hi  = 100*(exp(T.ci97_5)-1);
labs= erase(string(T.param),'beta:');

% 排序与截取
[~,ord] = sort(abs(m),'descend');
ord = ord(1:min(opt.Kshow,numel(ord)));
m=m(ord); lo=lo(ord); hi=hi(ord); labs=labs(ord);

% ——— 绘图 ———
fh=figure('Color','w','Units','pixels','Position',[80 80 540 680]);
ax=axes(fh); hold(ax,'on');

y = 1:numel(m);                                        % ❶ 递增向量
for i=1:numel(y)
    plot(ax,[lo(i) hi(i)],[y(i) y(i)],'-','Color',[0.15 0.15 0.15 0.9],'LineWidth',1.0);
end
plot(ax,m,y,'o','MarkerFaceColor',opt.Color,'MarkerEdgeColor','k','MarkerSize',5);

xline(ax,0,'-','Color',[0.6 0.6 0.6],'LineWidth',1);

set(ax,'YTick',y,'YTickLabel',labs,'YDir','reverse', ...   % ❷ 反转 y 轴方向
    'Box','off','TickDir','out','LineWidth',0.8,'FontName',opt.FontName,'FontSize',14);

xlabel(ax,'Relative change (%)','FontName',opt.FontName,'FontSize',13);
title(ax,panel_title,'FontWeight','normal','FontName',opt.FontName,'FontSize',14);

grid(ax,'on'); ax.GridColor=[0.9 0.9 0.9]; ax.GridAlpha=0.8;

% —— 自适应边界（鲁棒版） ——
% 仅用有限数计算范围；必要时回退到点估计 m
finiteMask = isfinite(lo) & isfinite(hi);
cand = [lo(finiteMask); hi(finiteMask)];
if isempty(cand)
    cand = m(isfinite(m));
end
% 如果依然为空，给一个保守兜底区间
if isempty(cand)
    xlo = -10; xhi = 10;
else
    xlo = min(cand); xhi = max(cand);
    if ~isfinite(xlo) || ~isfinite(xhi)
        xlo = -10; xhi = 10;
    end
end
% 若上下界相等或顺序不增，轻微扩一丢丢
if ~(xhi > xlo)
    epsPad = max(1, 0.1*max(1,abs(xhi)));
    xlo = xlo - epsPad;
    xhi = xhi + epsPad;
end
% 再加 5% padding
pad = 0.05*(xhi - xlo);
xlo = xlo - pad; xhi = xhi + pad;
% 最终保证严格递增
if xhi <= xlo, xhi = xlo + 1; end
xlim(ax, [xlo, xhi]);


% 导出
if ~isempty(opt.Out)
    outdir=fileparts(char(opt.Out)); if ~isfolder(outdir), mkdir(outdir); end
    exportgraphics(fh,[char(opt.Out) '.png'],'Resolution',opt.DPI);
    print(fh,[char(opt.Out) '.pdf'],'-dpdf','-painters');
end
end
