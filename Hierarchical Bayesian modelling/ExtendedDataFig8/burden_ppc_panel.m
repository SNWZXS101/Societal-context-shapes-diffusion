function burden_ppc_panel(results_root, target_dir, panel_title, varargin)
% BURDEN_PPC_PANEL (robust)
% PPC 对数散点：x=观测, y=预测中位数；误差条=10–90%预测区间；颜色=变异株
% 若缺失 data.json → 自动降级为单色绘制，不分变异株。
%
% 用法：
%   burden_ppc_panel('results_hbm','I_max', 'Posterior predictive: I_{max}', ...
%       'Colors',[0.82 0.33 0.30; 0.95 0.66 0.12; 0.12 0.47 0.71], ...
%       'Out','results_hbm/figs_burden/PanelC_Imax');

p=inputParser;
addParameter(p,'Colors',[0.82 0.33 0.30; 0.95 0.66 0.12; 0.12 0.47 0.71],@(x)isnumeric(x)&&size(x,2)==3);
addParameter(p,'Out','',@(s)ischar(s)||isstring(s));
addParameter(p,'FontName','Arial',@(s)ischar(s)||isstring(s));
addParameter(p,'DPI',400,@(x)isnumeric(x)&&x>0);
parse(p,varargin{:});
opt=p.Results;

ppc_path = fullfile(results_root,'burden',target_dir,'ppc.csv');
if ~isfile(ppc_path)
    error('未找到 %s', ppc_path);
end
ppc = readtable(ppc_path);

% ------- 读取 variant_id（若 data.json 存在） -------
dj_path = fullfile(results_root,'burden',target_dir,'data.json');
has_json = isfile(dj_path);
if has_json
    try
        dj = jsondecode(fileread(dj_path));
    catch
        warning('[PPC] data.json 读取失败，将不区分变异株着色。');
        has_json = false;
    end
end

x = double(ppc.y_obs(:));
y = double(ppc.yhat_p50(:));
hasBands = all(ismember({'yhat_p10','yhat_p90'}, ppc.Properties.VariableNames));
if hasBands
    ylo = double(ppc.yhat_p10(:));
    yhi = double(ppc.yhat_p90(:));
end

% ------- 清理非法/缺失 -------
ok = isfinite(x) & isfinite(y);
if hasBands
    ok = ok & isfinite(ylo) & isfinite(yhi);
end
x = x(ok); y = y(ok);
if hasBands, ylo = ylo(ok); yhi = yhi(ok); end

% ------- 变异株着色 -------
if has_json && isfield(dj,'variant_id') && numel(dj.variant_id)>=numel(ok)
    var_id_full = double(dj.variant_id(:));
    var_id = var_id_full(ok);
    use_groups = true;
else
    % 找不到 data.json 或无 variant_id → 降级为单色
    var_id = ones(numel(x),1);
    use_groups = false;
    if ~has_json
        warning('[PPC] %s 缺失，已降级为单色绘制；请确认目标子目录名称是否正确（如 I_max vs I_max_{obs}）。', dj_path);
    end
end
V = max(var_id);
cols = opt.Colors;
if size(cols,1) < V
    % 若颜色不够，补充重复或生成简易调色
    cols = lines(V);
end

% ------- 画图 -------
% ------- 画图 -------
fh=figure('Color','w','Units','pixels','Position',[80 80 520 480]);
ax=axes(fh); hold(ax,'on');

hscat = gobjects(V,1);                           % <— 用来记录每个变异株的散点句柄
for v=1:V
    idx = (var_id==v);
    if ~any(idx), continue; end

    if hasBands
        % 竖向误差条（不进图例）
        for i=find(idx').'
            hbar = plot(ax, [x(i) x(i)], [ylo(i) yhi(i)], '-', ...
                'Color',[0 0 0 0.15], 'LineWidth',1);
            % 这行把误差条从图例中隐藏
            for hb = hbar(:)'
                hb.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end
        end
    end

    % 散点（进图例）
    hscat(v) = scatter(ax, x(idx), y(idx), 28, cols(v,:), 'filled', ...
        'MarkerEdgeColor','k','MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.9, ...
        'DisplayName', sprintf('Variant %d', v));
end


% 45°线 + 轴
posx = x(x>0); posy = y(y>0);
if isempty(posx), posx = x; end
if isempty(posy), posy = y; end
epsx = max(min(posx), min([posx;posy])*0.8);
epsx = max(epsx, 1e-6);                 % log 轴兜底
mx   = max([x; y; hasBands*max(yhi)]);
mx   = mx * 1.15;

plot(ax,[epsx mx],[epsx mx],'--','Color',[0.3 0.3 0.3],'LineWidth',1);
set(ax,'XScale','log','YScale','log','Box','off','TickDir','out','LineWidth',0.8,...
    'FontName',opt.FontName,'FontSize',13);
grid(ax,'on'); ax.GridColor=[0.9 0.9 0.9]; ax.GridAlpha=0.85;

% 轴范围（严格递增兜底）
loBound = min([x; y; hasBands*ylo]); if ~isfinite(loBound) || loBound<=0, loBound = epsx; end
hiBound = max([x; y; hasBands*yhi]); if ~isfinite(hiBound) || hiBound<=loBound, hiBound = loBound*10; end
xlim(ax,[loBound, hiBound]); ylim(ax,[loBound, hiBound]);

xlabel(ax,'Observed','FontName',opt.FontName,'FontSize',13);
ylabel(ax,'Predicted median','FontName',opt.FontName,'FontSize',13);
title(ax,panel_title,'FontWeight','normal','FontName',opt.FontName,'FontSize',14);

% 角标：Spearman & 覆盖率
rho = corr(log(x+eps), log(y+eps), 'type','Spearman','rows','pairwise');
if hasBands
    cov10_90 = mean( (x>=ylo) & (x<=yhi) );
else
    cov10_90 = NaN;
end
text(ax,0.02,0.98, sprintf('\\rho=%.3f   cov_{10–90}=%.3f', rho, cov10_90), ...
    'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top', ...
    'FontSize',13,'Color',[0.1 0.1 0.1]);

% 图例（只有在有多组时显示）
if use_groups && V>1
    hleg = hscat(isgraphics(hscat));
    if ~isempty(hleg)
        lgd = legend(hleg, 'Location','southoutside','Orientation','horizontal');
        lgd.Box='off'; lgd.ItemTokenSize=[10 10];
    end
end


% 导出
if ~isempty(opt.Out)
    outdir=fileparts(char(opt.Out)); if ~isfolder(outdir), mkdir(outdir); end
    exportgraphics(fh,[char(opt.Out) '.png'],'Resolution',opt.DPI);
    print(fh,[char(opt.Out) '.pdf'],'-dpdf','-painters');
end
end
