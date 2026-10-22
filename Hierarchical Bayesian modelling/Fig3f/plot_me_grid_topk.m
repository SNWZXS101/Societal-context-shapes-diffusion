function h = plot_me_grid_topk(results_root, varargin)
% PLOT_ME_GRID_TOPK
% Fig.2b — Grid: columns = targets; rows = per-column Top-K covariates.
% Each tile shows 10–90% band + mean line (inverse-transformed response).
%
% Requirements:
%  results_root/M0_main/<target>/marginal_effects.csv with:
%    covariate, grid_std, yhat_mean, yhat_p10, yhat_p50, yhat_p90

% -------- params --------
p = inputParser;
addParameter(p,'Targets',{},@(c)iscell(c)||isstring(c));
addParameter(p,'TopK',6,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'SelectBy','amplitude',@(s)ischar(s)||isstring(s));
addParameter(p,'TargetPalette','maina8',@(s)ischar(s)||isstring(s));
addParameter(p,'FontSize',9,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'BandAlpha',0.20,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
addParameter(p,'MeanLineWidth',1.1,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'Title','',@(s)ischar(s)||isstring(s));
addParameter(p,'XRange',[-2 2],@(v)isnumeric(v)&&numel(v)==2&&v(1)<v(2));
addParameter(p,'ShowYTicks',true,@islogical);          % keep control of y-ticks
addParameter(p,'SameYPerColumn',false,@islogical);     % optional unify per-column y-limits

addParameter(p,'YLabelWrap',true,@islogical);          % 是否换行
addParameter(p,'YLabelMaxChars',12,@(x)isnumeric(x)&&isscalar(x)&&x>=4); % 每行最多字符数
addParameter(p,'YLabelRotation',90,@(x)isnumeric(x)&&isscalar(x));       % 旋转角度（默认竖排）


parse(p,varargin{:});
A = p.Results;

% -------- targets & colours --------
S = load(fullfile(results_root,'cache','hbm_cache.mat')); cache=S.cache;
targets = A.Targets; if isempty(targets), targets = cellstr(string(cache.active_targets)); end
nT = numel(targets);
CT = smf_colorlib(A.TargetPalette,'K',nT);

% -------- read marginal effects --------
ME = struct();
for ti=1:nT
    tar = string(targets{ti});
    fcsv = fullfile(results_root,'M0_main',char(tar),'marginal_effects.csv');
    if ~isfile(fcsv), error('Missing %s', fcsv); end
    T = readtable(fcsv);
    need = {'covariate','grid_std','yhat_mean','yhat_p10','yhat_p50','yhat_p90'};
    assert(all(ismember(need, string(T.Properties.VariableNames))), ...
        'marginal_effects.csv missing required columns for %s', tar);
    T = T(:, need);
    T.Properties.VariableNames = {'covariate','grid_std','mean','p10','p50','p90'};
    ME.(matlab.lang.makeValidName(char(tar))) = T;
end

% -------- select Top-K per column (target) --------
RowCovars = cell(nT,1);
for col=1:nT
    tar = string(targets{col});
    T = ME.(matlab.lang.makeValidName(char(tar)));
    covs = unique(T.covariate,'stable');
    sc = zeros(numel(covs),1);
    for i=1:numel(covs)
        Ci = sortrows(T(strcmp(T.covariate,covs{i}),:),'grid_std');
        switch lower(string(A.SelectBy))
            case 'amplitude',      sc(i) = max(Ci.mean)-min(Ci.mean);
            case 'range_p90p10',   sc(i) = max(Ci.p90)-min(Ci.p10);
            case 'l2',             sc(i) = sqrt(mean(Ci.mean.^2));
            case 'absmean',        sc(i) = mean(abs(Ci.mean));
            case 'slope',          sc(i) = abs(Ci.mean(end)-Ci.mean(1));
            otherwise, error('Unknown SelectBy: %s', A.SelectBy);
        end
    end
    [~,ord] = sort(sc,'descend');
    k = min(A.TopK,numel(covs));
    RowCovars{col} = cellstr(covs(ord(1:k)));
end

% optional per-column y-limits
YlimCol = cell(1, nT);
if A.SameYPerColumn
    for col = 1:nT
        tar = string(targets{col});
        T   = ME.(matlab.lang.makeValidName(char(tar)));
        covs_col = RowCovars{col};
        ymin = inf; ymax = -inf;
        for r=1:numel(covs_col)
            C = T(strcmp(T.covariate, covs_col{r}),:);
            if ~isempty(C)
                ymin = min([ymin; C.p10(:)]);
                ymax = max([ymax; C.p90(:)]);
            end
        end
        if isfinite(ymin) && isfinite(ymax) && ymin<ymax
            pad = 0.04*(ymax - ymin + eps);
            YlimCol{col} = [ymin-pad, ymax+pad];
        else
            YlimCol{col} = [];
        end
    end
end

% -------- layout & draw --------
fig = figure('Color','w','Name','Grid (TopK per column)');
tlo = tiledlayout(fig, A.TopK, nT, 'Padding','compact','TileSpacing','compact');
ax_mat = gobjects(A.TopK, nT);

for r=1:A.TopK
    for col=1:nT
        ax = nexttile(tlo, (r-1)*nT + col);
        hold(ax,'on'); box(ax,'off');
        set(ax,'FontName','Helvetica','FontSize',A.FontSize, ...
            'LineWidth',0.75,'TickDir','out');
        grid(ax,'on'); ax.GridColor=[0.95 0.95 0.95]; ax.GridAlpha=1.0;

        tar  = string(targets{col});
        colC = CT(col,:);
        T    = ME.(matlab.lang.makeValidName(char(tar)));
        covs_col = RowCovars{col};

        haveData = false; x = [];
        if r <= numel(covs_col)
            cname = string(covs_col{r});
            C = sortrows(T(strcmp(T.covariate, cname), :), 'grid_std');
            if ~isempty(C)
                x = C.grid_std(:);
                Xp = [x; flipud(x)];
                Yp = [C.p10(:); flipud(C.p90(:))];

                patch('XData',Xp,'YData',Yp, ...
                    'FaceColor',colC,'FaceAlpha',A.BandAlpha, ...
                    'EdgeColor','none','Parent',ax);
                plot(ax, x, C.mean(:), '-', 'Color',[1 1 1], ...
                    'LineWidth', A.MeanLineWidth);

                haveData = true;

                % ★ 每个小面板的纵轴标签 = 该面板的协变量名
                ylab_str = format_covariate_label(cname, A);  % ★ 新：自动换行 + TeX 转义
                ylabel(ax, ylab_str, 'Interpreter','tex', ...
                    'FontSize', A.FontSize+0.5, 'FontWeight','bold', ...
                    'Rotation', A.YLabelRotation, ...
                    'HorizontalAlignment','center');

            end
        end

        % top row: column title = target name (colour-matched)
        if r==1
            title(ax, format_target_label(tar), 'Interpreter','tex', ...
                'Color', colC, 'FontWeight','bold', 'FontSize', A.FontSize+1);
        end

        % y-ticks on/off
        if ~A.ShowYTicks, yticks(ax, []); end

        % per-column unified y-limits (optional)
        if A.SameYPerColumn && ~isempty(YlimCol{col})
            ylim(ax, YlimCol{col});
        end

        % x-limits = union of data and XRange
        if haveData
            xmin = min([x(:); A.XRange(1)]);
            xmax = max([x(:); A.XRange(2)]);
        else
            xmin = A.XRange(1); xmax = A.XRange(2);
        end
        if ~isfinite(xmin) || ~isfinite(xmax) || xmin>=xmax
            xmin = A.XRange(1); xmax = A.XRange(2);
        end
        xlim(ax, [xmin xmax]);

        % bottom row: x label
        if r==A.TopK
            xlabel(ax, 'z (standardized covariate)', 'Interpreter','tex');
        end

        ax_mat(r,col) = ax;
    end
end

% super title
if strlength(string(A.Title))>0
    sgtitle(tlo, string(A.Title), 'Interpreter','tex', ...
        'FontSize', A.FontSize+5, 'FontWeight','bold');
end

h = struct('fig',fig,'ax',ax_mat,'targets',{targets},'selected',RowCovars);
end

% ---- helper ----
function s = format_target_label(t)
switch char(t)
    case 'delta_early',      s = '\Delta_{early}';
    case 'kappa_cross',      s = '\kappa_{cross}';
    case 'rho_eps',          s = '\rho_{\epsilon}';
    case 'tf1_mean_early',   s = 'f_{1,early}';
    case 'tf2_peak',         s = 'f_{2,peak}';
    case 't_tf2_peak',       s = 't_{f2peak}';
    otherwise,               s = char(t);
end
end

function out = format_covariate_label(cname, A)
% 自动换行 + TeX 转义。优先在空格/下划线处分行。
s = string(cname);

% 让下划线也能作为“词”的边界，便于在它附近换行
s = replace(s, "_", " _ ");

% 分词（按任意空白）
tokens = split(s);                      % string 向量
tokens = tokens(strlength(tokens)>0);   % 去掉空 token

% 逐个拼接，超过最大长度就换新行
if A.YLabelWrap
    lines = strings(0,1);
    cur = "";
    for i = 1:numel(tokens)
        t = tokens(i);
        if strlength(cur) == 0
            cand = t;
        else
            cand = cur + " " + t;
        end

        if strlength(cand) <= A.YLabelMaxChars
            cur = cand;
        else
            % 先收当前行
            if strlength(cur) > 0
                lines(end+1,1) = cur; %#ok<AGROW>
            end
            % 若单词本身超长，硬切分
            if strlength(t) > A.YLabelMaxChars
                tt = char(t); k = A.YLabelMaxChars;
                while ~isempty(tt)
                    piece = string(tt(1:min(k,end)));
                    lines(end+1,1) = piece; %#ok<AGROW>
                    tt = tt(min(k,end)+1:end);
                end
                cur = "";
            else
                cur = t;
            end
        end
    end
    if strlength(cur) > 0
        lines(end+1,1) = cur;
    end
    s = strjoin(lines, newline);        % 用换行连接
else
    s = strjoin(tokens, " ");
end

% —— TeX 转义（仅标题/标签用到 'tex' 时需要；这里始终转义更稳）——
s = replace(s, "\", "\\");
s = replace(s, "_", "\_");
s = replace(s, "%", "\%");
s = replace(s, "^", "\^");
s = replace(s, "{", "\{");
s = replace(s, "}", "\}");

out = s;
end
