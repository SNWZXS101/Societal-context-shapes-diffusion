function create_fig_srs_ranking(scores_root, varargin)
% CREATE_FIG_SRS_RANKING  30省 SRS 横向排名图（同色系序列色，调用 colormap_factory）
% 用法：
%   create_fig_srs_ranking('scores', ...
%       'style','blue', ...                 % 'orange'|'purple'|'blue'|'green'|'red'|'yellow'|'teal'
%       'export_pdf',true,'out_pdf','Fig_SRS_Ranking.pdf');

p = inputParser; p.FunctionName = mfilename;
addRequired(p,'scores_root',@(s)ischar(s)||isstring(s));
addParameter(p,'style','blue',@(s)ischar(s)||isstring(s));
addParameter(p,'export_pdf',false,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'out_pdf','Fig_SRS_Ranking.pdf',@(s)ischar(s)||isstring(s));
parse(p,scores_root,varargin{:}); A=p.Results;

% ---- 读 SRS ----
T = readtable(fullfile(A.scores_root,'map_ready_SRS.csv'));  % province, RRI_VA, MEI_VA, SRS, rank
T.province = string(T.province);

% ---- 排序（SRS 降序）----
[~,ord] = sort(T.SRS, 'descend','MissingPlacement','last');
T = T(ord,:);

% ---- 颜色：调用你的 colormap_factory(style, N)；Top=最深 ----
N = height(T);
Cseq = colormap_factory(A.style, N);  % 从浅到深
barColor = flipud(Cseq);              % 让 rank1 用最深

% ---- 画布 ----
fig = figure('Color','w','Position',[120 100 900 900]);
movegui(fig,'center'); plotedit(fig,'on');

ax  = axes('Parent',fig,'Position',[0.22 0.08 0.58 0.88],'Color','w');
hold(ax,'on'); box(ax,'on'); grid(ax,'on');
set(ax,'GridAlpha',0.11,'GridColor',[0 0 0],'MinorGridAlpha',0.05);

% ---- 横向条形 ----
y = 1:N;
bh = barh(ax, y, T.SRS, 'FaceColor','flat', 'BarWidth',0.72, ...
          'EdgeColor',[0.12 0.12 0.12], 'LineWidth',0.6);
bh.CData = barColor;

% 纵轴用完整省名；从上到下排名递增（Top 在最上）
set(ax,'YLim',[0.5 N+0.5],'YTick',y,'YTickLabel',T.province,'YDir','reverse');
xlabel(ax, 'SRS (0–100)');
title(ax, 'SRS Ranking of 30 Provinces');

% ---- 只在柱子尾部标注数值（不重复省名）----
xmax = max(T.SRS,[],'omitnan'); xlim(ax,[0 max(100, xmax*1.06)]);
for i=1:N
    text(ax, T.SRS(i)+0.8, y(i), sprintf('%.1f', T.SRS(i)), ...
        'HorizontalAlignment','left','VerticalAlignment','middle', ...
        'FontSize',12,'Color',[0.15 0.15 0.15], 'Interpreter','none');
end


% % —— 紧凑的“排名色带”示意（深=高SRS，浅=低SRS）——
% cax = axes('Parent',fig,'Position',[0.83 0.10 0.015 0.82], ...
%            'Box','off','XTick',[],'YTick',[],'Color','none');
% 
% Cstrip = flipud(Cseq);                       % 关键：翻转→顶端为深色
% imagesc(cax, permute(Cstrip,[1 3 2]));       % 竖向色带：上深下浅
% set(cax,'YDir','normal');                    % 保持顶端为深
% 
% axis(cax,'off');
% text(cax, 1.8, 0.95,'High','Units','normalized','Rotation',90, ...
%      'HorizontalAlignment','center','FontSize',8,'Color',[0.25 0.25 0.25]);
% text(cax, 1.8, 0.05, 'Low','Units','normalized','Rotation',90, ...
%      'HorizontalAlignment','center','FontSize',8,'Color',[0.25 0.25 0.25]);
% title(cax, 'Rank color', 'FontSize',14,'Color',[0.25 0.25 0.25]);


% ---- 导出 ----
% if A.export_pdf
%     set(fig,'Renderer','painters');
%     print(fig, A.out_pdf, '-dpdf','-bestfit');
%     fprintf('[export] saved: %s\n', A.out_pdf);
% end
end
