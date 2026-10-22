function h = plot_volcano_overlay_labels(results_root, varargin)
% PLOT_VOLCANO_OVERLAY_LABELS
% 单图叠加火山图（6 目标），每目标自动标注“最大正/最大负”两点（共12个标签）
%
% 示例：
%   h = plot_volcano_overlay_labels('results_hbm', ...
%        'Targets', {'delta_early','kappa_cross','rho_eps','tf1_mean_early','tf2_peak','t_tf2_peak'}, ...
%        'UseDraws', true, ...
%        'TargetPalette','maina8', ...
%        'FontSize', 9, ...            % 轴与标题字号
%        'MarkerSize', 64, ...         % 所有散点大小
%        'LabelFontSize', 8, ...       % 标签文本框字号
%        'LeaderLineWidth', 1.0, ...   % 指向线粗细
%        'Title', 'Fig.2a  Volcano (All targets overlaid)', ...
%        'YMax', 3.2, 'XThresh', [-0.2 0.2]);  % 竖向阈值线可选
%
% 参数（Name–Value）：
%   'Targets'        : 目标列表（默认从 cache.active_targets 读取）
%   'UseDraws'       : 是否用 draws_raw 精确计算 P_same（默认 true）
%   'TargetPalette'  : 目标配色（默认 'maina8' 取前6色）
%   'FontSize'       : 轴刻度/标题/标签的基准字号（默认 8）
%   'MarkerSize'     : 所有散点大小（默认 48）
%   'LabelFontSize'  : 标签文本框字号（默认 8）
%   'LeaderLineWidth': 标签指向线宽（默认 0.8）
%   'YMax'           : 纵轴上限（默认 3.2；0.95/0.99/0.999≈1.301/2/3）
%   'XThresh'        : 竖向阈值线数组（例如 [-0.2 0.2]，默认 [] 不画）
%   'Title'          : 图标题（tex）
%
% 依赖：
%   results_root/cache/hbm_cache.mat
%   results_root/M0_main/<target>/summary.csv
%   results_root/M0_main/<target>/fe_meta.mat
%   （可选）results_root/M0_main/<target>/draws_raw.mat
%   smf_colorlib.m

% ---------------- 参数 ----------------
p = inputParser;
addParameter(p,'Targets',{},@(c)iscell(c)||isstring(c));
addParameter(p,'UseDraws',true,@islogical);
addParameter(p,'TargetPalette','maina8',@(s)ischar(s)||isstring(s));
addParameter(p,'FontSize',8,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'MarkerSize',48,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'LabelFontSize',8,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'LeaderLineWidth',0.8,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'YMax',3.2,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'XThresh',[],@(v)isnumeric(v)&&isvector(v));
addParameter(p,'Title','',@(s)ischar(s)||isstring(s));
addParameter(p,'LeaderMode','annotation',@(s)ischar(s)||isstring(s));
parse(p,varargin{:});
A = p.Results;

% ---------------- 数据准备 ----------------
S = load(fullfile(results_root,'cache','hbm_cache.mat')); cache = S.cache;
targets = A.Targets; if isempty(targets), targets = cellstr(string(cache.active_targets)); end
Ctar = smf_colorlib(A.TargetPalette, 'K', numel(targets));


% ===== NEW: 变量组映射 + per-target 容器 =====
% groupMap: fe_name -> group（从 cache 中尽可能提取）
groupMap = containers.Map;
try
    % 常见字段名的容错
    if isfield(cache,'fe_names_orig') && isfield(cache,'fe_groups')
        k = string(cache.fe_names_orig(:));
        v = string(cache.fe_groups(:));
        groupMap = containers.Map(k, v);
    elseif isfield(cache,'fe_names') && isfield(cache,'fe_groups')
        k = string(cache.fe_names(:));
        v = string(cache.fe_groups(:));
        groupMap = containers.Map(k, v);
    end
catch
end

% SumByTarget: 每个目标一个 table（后面用于统计/导出）
SumByTarget = containers.Map;



ALL = table();   % name, mean, sd, P_same, S, target
FE_LIST = cell(numel(targets),1);

for ti=1:numel(targets)
    tar = string(targets{ti});
    ydir = fullfile(results_root,'M0_main', char(tar));
    Tsum = readtable(fullfile(ydir,'summary.csv'));
    M = load(fullfile(ydir,'fe_meta.mat')); meta = M.meta;
    fe_names = string(meta.fe_names(:));
    FE_LIST{ti} = fe_names;

    % 只取 beta 行，并恢复名称
    isb = startsWith(string(Tsum.param),'beta','IgnoreCase',true);
    Tb = Tsum(isb,:);
    names = strings(height(Tb),1);
    for i=1:height(Tb)
        pi = string(Tb.param(i));
        if contains(pi,':')
            names(i) = extractAfter(pi,':');
        else
            tok = regexp(pi,'\d+$','match','once');
            if ~isempty(tok)
                k = str2double(tok);
                if k>=1 && k<=numel(fe_names), names(i)=fe_names(k); else, names(i)="beta_"+tok; end
            else
                names(i) = pi;
            end
        end
    end

    % 证据：P_same
    if A.UseDraws && isfile(fullfile(ydir,'draws_raw.mat'))
        P_same = psame_from_draws(fullfile(ydir,'draws_raw.mat'), fe_names, names);
    else
        z = abs(Tb.mean)./max(Tb.sd, eps);
        P_same = normcdf(z);
    end
    Sval = -log10(1 - min(max(P_same, eps), 1-1e-9));

    Tthis = table(names, Tb.mean, Tb.sd, P_same, Sval, repmat(tar,height(Tb),1), ...
        'VariableNames', {'name','mean','sd','P_same','S','target'});

    % ===== NEW: 附加变量组列，并放入容器 =====
    G = strings(height(Tthis),1);
    for ii = 1:height(Tthis)
        key = string(Tthis.name(ii));
        if ~isempty(groupMap) && isKey(groupMap, key)
            G(ii) = string(groupMap(key));
        else
            G(ii) = "Unknown";
        end
    end
    Tthis.group = G;

    % 放入 per-target 容器（用于后续统计导出）
    SumByTarget(char(tar)) = Tthis;


    ALL = [ALL; Tthis]; %#ok<AGROW>
end

% ---------------- 作图 ----------------
fig = figure('Color','w','Name','Fig2a Volcano (overlay)');
ax  = axes(fig); hold(ax,'on'); box(ax,'off');
set(ax,'FontName','Helvetica','FontSize',A.FontSize,'LineWidth',0.75,'TickDir','out');

% 竖向阈值线（可选）
if ~isempty(A.XThresh)
    for xk = A.XThresh(:).'
        plot(ax, [xk xk], ylim(ax), ':', 'Color',[0.7 0.7 0.7], 'LineWidth',0.75);
    end
end
% 横向证据线（0.95/0.99/0.999）
thr = [0.95 0.99 0.999]; Sthr = -log10(1-thr);
for k=1:numel(Sthr)
    yk = Sthr(k);
    plot(ax, xlim(ax), [yk yk], ':', 'Color',[0.7 0.7 0.7], 'LineWidth',0.75);
end

% 逐目标绘制散点（未标注：无描边）
legH = gobjects(numel(targets),1);
for ti=1:numel(targets)
    Tt = ALL(strcmp(string(ALL.target), string(targets{ti})), :);
    col = Ctar(ti,:);
    hsc = scatter(ax, Tt.mean, Tt.S, A.MarkerSize, ...
        'Marker','o', 'MarkerFaceColor', col, ...
        'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.95, ...
        'DisplayName', string(targets{ti}));
    legH(ti)=hsc;
end
legend(ax, legH, 'Location','eastoutside', 'Interpreter','tex', 'Box','off');

% 自动寻找并标注：每目标最大正/最大负各 1 个
for ti=1:numel(targets)
    tar = string(targets{ti});
    Tt = ALL(strcmp(string(ALL.target), tar), :);
    col = Ctar(ti,:);

    % 最大正 / 最大负（按 mean）
    [~, ipos] = max(Tt.mean);
    [~, ineg] = min(Tt.mean);
    pick = unique([ipos ineg]);

    for idx = pick
        x = Tt.mean(idx); y = Tt.S(idx);
        % 加“有描边”的点（覆盖掉无描边版本）
        scatter(ax, x, y, A.MarkerSize, 'o', ...
            'MarkerFaceColor', col, 'MarkerEdgeColor','k', ...
            'LineWidth', 0.8, 'MarkerFaceAlpha', 1.0, 'HandleVisibility','off');

        % 标签锚点（数据坐标）：向外偏移
        dx = 0.02 * max(1, max(abs(ALL.mean)));
        dy = 0.06;
        xlab = x + sign(x + (x==0))*dx;
        ylab = y + dy;

        if strcmpi(A.LeaderMode,'annotation')
            % ---- annotation 版本（可拖动端点/文本框）----
            % 把数据坐标转到 figure-normalized
            [x1n, y1n] = axdata2fignorm(ax, x,    y);
            [x2n, y2n] = axdata2fignorm(ax, xlab, ylab);

            % 引线（可拖动端点）
            annotation(ax.Parent, 'line', [x1n x2n], [y1n y2n], ...
                'Color', col, 'LineWidth', A.LeaderLineWidth,  'HandleVisibility','off');

            % 文本框（可拖动/缩放），位置用 [x y w h]，先给个很小的 w,h 再自适应
            tb = annotation(ax.Parent, 'textbox', [x2n y2n 0.001 0.001], ...
                'String', string(Tt.name(idx)), 'Interpreter','tex', ...
                'FitBoxToText','on', 'BackgroundColor', col, 'EdgeColor','k', ...
                'Color', text_contrast_color(col), 'FontSize', A.LabelFontSize, ...
                'VerticalAlignment','bottom', 'HorizontalAlignment','left',...
                'Margin', 3);
        else
            % ---- 退回 axes-objects 版本（不可拖端点，只可整体拖）----
            line(ax, [x, xlab], [y, ylab], 'Color', col, ...
                'LineWidth', A.LeaderLineWidth, 'Clipping','on', 'HandleVisibility','off');
            text(ax, xlab, ylab, string(Tt.name(idx)), ...
                'Interpreter','tex', 'FontSize', A.LabelFontSize, ...
                'FontWeight','bold', 'Color', text_contrast_color(col), ...
                'BackgroundColor', col, 'Margin', 3, 'EdgeColor','k', ...
                'HorizontalAlignment','left', 'VerticalAlignment','bottom', ...
                'Clipping','on');
        end
    end
end

% 轴与标题
xlabel(ax, '\beta (posterior mean, transformed scale)', 'Interpreter','tex', 'FontSize', A.FontSize);
ylabel(ax, 'Evidence  S = -log_{10}(1 - P_{same})', 'Interpreter','tex', 'FontSize', A.FontSize);
grid(ax,'on'); ax.GridColor=[0.92 0.92 0.92]; ax.GridAlpha=1.0;

ylim(ax, [0, A.YMax]);
mx = max(1e-6, max(abs(ALL.mean)));
xlim(ax, mx*[-1.08 1.08]);

if strlength(string(A.Title))>0
    title(ax, string(A.Title), 'Interpreter','tex', 'FontSize', A.FontSize+2, 'FontWeight','bold');
end


% ===== NEW (fixed): 统计导出 —— 每目标 Top +/- 与分组阈值计数 =====
% ===== Export a concise volcano report (top +/- per target) =====
try
    rows = {};  % cell rows we will convert to table
    thr_list = [0.95 0.99 0.999];

    for ti = 1:numel(targets)
        tar = string(targets{ti});
        Tt  = ALL(strcmp(string(ALL.target), tar), :);

        % counts of coefficients with P_same >= thresholds (for this target)
        cnt95  = sum(Tt.P_same >= thr_list(1));
        cnt99  = sum(Tt.P_same >= thr_list(2));
        cnt999 = sum(Tt.P_same >= thr_list(3));

        % locate top positive / top negative (by posterior mean)
        [~, ipos] = max(Tt.mean);
        [~, ineg] = min(Tt.mean);

        pick_idx  = [ipos ineg];
        pick_sign = {'pos','neg'};

        for k = 1:numel(pick_idx)
            ii = pick_idx(k);
            if isempty(ii) || isnan(ii), continue; end

            rows(end+1,1:9) = { ...
                char(tar), ...                           % target
                pick_sign{k}, ...                        % sign
                char(Tt.name(ii)), ...                   % covariate
                Tt.mean(ii), ...                         % beta_mean (transformed scale)
                Tt.P_same(ii), ...                       % P_same
                Tt.S(ii), ...                            % S = -log10(1-P_same)
                cnt95, cnt99, cnt999 ...                 % counts above thresholds
            };
        end
    end

    if ~isempty(rows)
        Report = cell2table(rows, 'VariableNames', { ...
            'target','sign','covariate','beta_mean','P_same','S', ...
            'n_Psame_ge_0p95','n_Psame_ge_0p99','n_Psame_ge_0p999'});

        out_dir  = fullfile(results_root,'M0_main');
        if exist(out_dir,'dir') ~= 7, out_dir = results_root; end
        out_csv  = fullfile(out_dir,'volcano_overlay_report.csv');
        writetable(Report, out_csv);

        % console print (brief)
        fprintf('[volcano] Report saved: %s\n', out_csv);
        disp(Report);

        % also return in handle
        try
            h.report = Report;  % h is your function output struct
        catch
            % ignore if h not yet defined
        end
    end
catch ME
    warning(sprintf('Volcano report export skipped: %s', ME.message));
end









% 输出句柄
h = struct('fig',fig,'ax',ax,'data',ALL);
end

% ---------- 工具：从 draws_raw 计算 P_same ----------
function P_same = psame_from_draws(draws_mat_path, fe_names_all, names_q)
L = load(draws_mat_path);
if ~isfield(L,'draws'), P_same = normcdf(zeros(numel(names_q),1)); return; end
T = L.draws; vars = string(T.Properties.VariableNames);
pat = "(?i)^beta[_:]*\s*(\d+)$";
idx_num = nan(numel(vars),1);
for k=1:numel(vars)
    tk = regexp(vars(k), pat, 'tokens','once');
    if ~isempty(tk), idx_num(k) = str2double(tk{1}); end
end
keep = ~isnan(idx_num);
if ~any(keep), P_same = normcdf(zeros(numel(names_q),1)); return; end
[~, ord] = sort(idx_num(keep)); cols = find(keep); cols = cols(ord);
B = table2array(T(:, cols));             % S_total × K
P_same = nan(numel(names_q),1);
for i=1:numel(names_q)
    nm = string(names_q(i)); k = find(fe_names_all==nm,1);
    if isempty(k) || k>size(B,2), continue; end
    bj = B(:,k); m = mean(bj,'omitnan'); sgn = sign(m); if sgn==0, sgn=1; end
    P_same(i) = mean(sign(bj)==sgn,'omitnan');
end
end

% ---------- 工具：给定底色返回对比度充足的文字颜色 ----------
function c = text_contrast_color(bg)
% 相对亮度阈值法（WCAG 近似）
L = 0.2126*bg(1) + 0.7152*bg(2) + 0.0722*bg(3);
if L > 0.6, c = [0 0 0]; else, c = [1 1 1]; end
end


function [xn, yn] = axdata2fignorm(ax, x, y)
% 把 axes 数据坐标 (x,y) 映射到 figure normalized 坐标 (xn, yn)
% 仅适用线性坐标轴（本图就是线性）
fig = ancestor(ax,'figure');
oldAU = ax.Units; ax.Units = 'normalized';
pos = ax.Position;  % [x0 y0 w h] in figure-normalized
ax.Units = oldAU;

xlimv = xlim(ax); ylimv = ylim(ax);
xn = pos(1) + (x - xlimv(1)) / (xlimv(2)-xlimv(1)) * pos(3);
yn = pos(2) + (y - ylimv(1)) / (ylimv(2)-ylimv(1)) * pos(4);
end


% —— 辅助：把 groupcounts 的 (group,count) 对齐到既定 level 顺序 —— 
% function v = align_counts(levels, gtbl)
%     v = zeros(numel(levels),1);
%     if isempty(gtbl), return; end
%     % gtbl 可能是 table(Group,Count) 或 (group,count)，统一大小写处理
%     vn = string(gtbl.Properties.VariableNames);
%     cName = vn(contains(vn,"Count",'IgnoreCase',true));
%     gName = vn(contains(vn,"Group",'IgnoreCase',true) | strcmpi(vn,"group"));
%     gc = table2array(gtbl(:, cName));
%     gg = string(table2array(gtbl(:, gName)));
%     for i = 1:numel(gg)
%         j = find(levels==gg(i),1);
%         if ~isempty(j), v(j) = gc(i); end
%     end
% end