function h = hbm_uq_plot_coeffs_grid(root_uq, target, ncols)
% HBM_UQ_PLOT_COEFFS_GRID
% Plot fixed-effect robustness (median & 95% CI) across specifications
% in ONE figure with multi-column layout. Only the FIRST column shows
% the y-axis ticks/labels (specification names); other columns hide them.
%
% Usage:
%   hbm_uq_plot_coeffs_grid('results_hbm_uq','delta_early',4)
%   kappa_cross   rho_eps  tf1_mean_early  tf2_peak  t_tf2_peak
%
% Inputs:
%   root_uq : UQ root dir (default 'results_hbm_uq')
%   target  : target name (default 'delta_early')
%   ncols   : number of columns (default 4)
%
% Requires files:
%   results_hbm_uq/<spec>/M0_main/<target>/summary.csv

if nargin<1 || isempty(root_uq), root_uq = 'results_hbm_uq'; end
if nargin<2 || isempty(target),  target  = 'delta_early';    end
if nargin<3 || isempty(ncols),   ncols   = 4;                 end

% Global render: LaTeX

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% Desired spec order (use only those that exist)
want_specs = {'base','no_HS','no_RS','no_RS_HS','priorsx0p5','priorsx2'};
spec_dirs  = {};
for i=1:numel(want_specs)
    if isfolder(fullfile(root_uq, want_specs{i}))
        spec_dirs{end+1} = want_specs{i}; %#ok<AGROW>
    end
end
if isempty(spec_dirs)
    error('No spec folders found in %s.', root_uq);
end

% -------- read summary.csv for each spec and collect betas --------
all_rows = {};
for s = 1:numel(spec_dirs)
    sumcsv = fullfile(root_uq, spec_dirs{s}, 'M0_main', target, 'summary.csv');
    if exist(sumcsv,'file') ~= 2
        warning('Missing: %s (skip)', sumcsv); continue;
    end
    T = readtable(sumcsv);
    if ~ismember('param', T.Properties.VariableNames), continue; end

    pnames = local_cellstr(T.param);
    keep   = startsWith(pnames, 'beta:');
    if ~any(keep), continue; end
    T = T(keep, :);

    needCols = {'median','ci2_5','ci97_5'};
    for c = 1:numel(needCols)
        if ~ismember(needCols{c}, T.Properties.VariableNames)
            error('%s missing column %s', sumcsv, needCols{c});
        end
    end

    bnames   = regexprep(pnames(keep), '^beta:\s*', '');
    S = table();
    S.spec    = repmat({spec_dirs{s}}, height(T), 1);
    S.spec_lx = repmat({lx_spec(spec_dirs{s})}, height(T), 1);
    S.coef    = bnames(:);
    S.coef_lx = cellfun(@lx_coefname, bnames(:), 'uni', 0);
    S.med     = T.median(:);
    S.lo      = T.ci2_5(:);
    S.hi      = T.ci97_5(:);
    all_rows{end+1} = S; %#ok<AGROW>
end
if isempty(all_rows), error('No beta rows found to plot.'); end
TAB = vertcat(all_rows{:});
TAB.spec    = local_cellstr(TAB.spec);
TAB.spec_lx = local_cellstr(TAB.spec_lx);
TAB.coef    = local_cellstr(TAB.coef);
TAB.coef_lx = local_cellstr(TAB.coef_lx);

% Unique coeffs (row facets)
coefs   = unique(TAB.coef, 'stable');
coefs_l = local_pick_first_per_key(TAB.coef, TAB.coef_lx);
nC = numel(coefs);

% Grid layout
ncols = max(1, round(ncols));
nrows = ceil(nC / ncols);

% Create figure/tiledlayout
h = figure('Color','w','Position',[60 60 1400 max(500, 230*nrows)]);
tl = tiledlayout(h, nrows, ncols, 'TileSpacing','compact','Padding','compact');

% 总标题（target 是字符串，如 'delta_early'）
tlab = lx_target_single(target);  % 返回片段，不含 $$
% tstr = sprintf(['Fixed effects robustness across specifications ' ...
%                 '(target: $\\mathtt{%s}$, i.e., $%s$)'], ...
%                 strrep(target,'_','\\_'), tlab);
% title(tl, tstr, 'Interpreter','latex', 'FontSize', 14);

title(tl, ['Fixed effects robustness across specifications (target: $' tlab '$)'], ...
      'Interpreter','latex', 'FontSize', 18);


% For aligning specs across panels, pre-build spec labels/positions
spec_lx_ordered = cellfun(@lx_spec, spec_dirs, 'uni', 0);

for k = 1:nC
    r = ceil(k / ncols);
    c = k - (r-1)*ncols;

    ax = nexttile(tl, (r-1)*ncols + c);

    % data for this coef across specs
    rows = strcmp(TAB.coef, coefs{k});
    S = TAB(rows, :);
    % reorder by spec_dirs order
    [~, ord] = ismember(S.spec, spec_dirs);
    [~, o2]  = sort(ord);
    S = S(o2, :);

    y = 1:height(S);

    hold(ax, 'on');
    % horizontal intervals
    for j = 1:height(S)
        line(ax, [S.lo(j) S.hi(j)], [y(j) y(j)], 'LineWidth', 2, 'Color', [0.2 0.2 0.8]);
        plot(ax, S.med(j), y(j), 'o', 'MarkerFaceColor',[0.1 0.4 0.9], ...
             'MarkerEdgeColor','w', 'MarkerSize', 6);
    end
    xline(ax, 0, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
    grid(ax,'on'); ax.XGrid='on'; ax.YGrid='on'; ax.GridAlpha=0.12;

    % y-axis: only first column shows tick labels
    if c == 1
        set(ax, 'YDir','reverse', 'YTick', y, 'YTickLabel', spec_lx_ordered, ...
            'FontSize', 13);
        ylabel(ax, 'Specifications','FontSize',13);
    else
        set(ax, 'YDir','reverse', 'YTick', y, 'YTickLabel', [], 'FontSize', 13);
        ax.YColor = [0.85 0.85 0.85];
    end

    % x-axis label只给底部一行
    if r == nrows
        xlabel(ax, 'Posterior (median, 95\% CI)', 'FontSize', 14);
    end

    title(ax, coefs_l{k}, 'FontSize', 15, 'Interpreter','latex');
    hold(ax, 'off');
end
end

% ----------------- helpers -----------------
function s = lx_spec(name_raw)
switch char(name_raw)
    case 'base',        s = 'base';
    case 'no_HS',       s = 'no $HS$';
    case 'no_RS',       s = 'no $RS$';
    case 'no_RS_HS',    s = 'no $RS{+}HS$';
    case 'priorsx0p5',  s = '$\mathrm{priors}\!\times\!0.5$';
    case 'priorsx2',    s = '$\mathrm{priors}\!\times\!2$';
    otherwise,          s = strrep(char(name_raw),'_','\_');
end
end

function s = lx_coefname(raw)
t = char(raw);
% gentle cleanups for readability
t = strrep(t, '_{pc}', '\_{pc}');
t = strrep(t, '_per_1k', '\_per\_1k');
t = strrep(t, 'urban population density', 'urban density');
t = strrep(t, 'public transit volume', 'public transit');
t = strrep(t, 'hospital beds per 1k', 'beds per 1k');
t = regexprep(t, '_', '\_'); % escape remaining underscores
s = t;
end

function c = local_cellstr(x)
if iscellstr(x)
    c = x;
elseif iscell(x)
    c = cellfun(@local_tochar, x, 'uni', 0);
elseif iscategorical(x)
    c = cellstr(x);
elseif ischar(x)
    c = {x};
else
    try, c = cellstr(x); catch, c = {local_tochar(x)}; end
end
end

function s = local_tochar(y)
if ischar(y), s = y;
elseif isnumeric(y), s = num2str(y);
else, try, s = char(y); catch, s=''; end
end
end

function out = local_pick_first_per_key(keys, labels)
keys  = local_cellstr(keys);
labels= local_cellstr(labels);
[u, ia] = unique(keys, 'stable');
out = labels(ia);
for i=1:numel(out)
    if isempty(out{i}), out{i} = strrep(u{i},'_','\_'); end
end
end


% ===== 目标名 → LaTeX 标签（单个）=====
function lab = lx_target_single(key)
% 返回 *不含* 外层 $...$ 的 LaTeX 片段，便于在外层统一包裹
switch string(key)
    case "delta_early",     lab = '\Delta_{\mathrm{early}}';
    case "kappa_cross",     lab = '\kappa_{\mathrm{cross}}';
    case "rho_eps",         lab = '\rho_{\epsilon}';
    case "tf1_mean_early",  lab = 'f_{1}^{\mathrm{early}}';
    case "tf2_peak",        lab = 'f_{2}^{\mathrm{peak}}';
    case "t_tf2_peak",      lab = 't_{\mathrm{f2peak}}';
    otherwise,              lab = strrep(key,'_','\_');  % 兜底：转义下划线
end
end

% （可选）批量版本，和你给的一样，便于别处复用
function L = lx_targets(keys)
L = cell(size(keys));
for i = 1:numel(keys)
    L{i} = ['$' lx_target_single(keys{i}) '$'];
end
end
