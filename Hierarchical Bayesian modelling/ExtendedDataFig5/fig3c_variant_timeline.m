function fig3c_variant_timeline(varargin)
% FIG3C_VARIANT_TIMELINE
% Variant timeline bands + wave scatters （Fig.3c 高颜值实现 · 修正版）
% -------------------------------------------------------------------------
% 依赖：
%   - results_hbm/cache/hbm_cache.mat
%   - results_hbm/M0_main/<target>/ranef_variant.csv
%   - features_aug.mat 或 features.mat
%
% 样式：
%   - 三条时间带（Original/Delta/Omicron）：95% CI 半透明带 + 中位线
%   - 波散点：颜色=变异株，点径=规模（cases_total/duration_days/w_tail_14 ...）
%   - 六宫格：默认 2×3，对应六个 target
%
% 调用示例：
%   fig3c_variant_timeline('results_root','results_hbm', ...
%       'features_mat','features_aug.mat', ...
%       'targets', {'delta_early','kappa_cross','rho_eps','tf1_mean_early','tf2_peak','t_tf2_peak'}, ...
%       'point_size_metric','w_tail', ...        % 'cases_total' | 'duration_days' | 'w_tail'
%       'smooth_bands',false, ...
%       'export_pdf','Fig3c_variant_timeline.pdf');

% -------------------- 参数 --------------------
p = inputParser;
addParameter(p,'results_root','results_hbm',@ischar);
addParameter(p,'features_mat','features_aug.mat',@ischar);
addParameter(p,'targets',{},@(c)iscellstr(c)||isstring(c));
addParameter(p,'point_size_metric','cases_total',@(s)ischar(s)||isstring(s));
addParameter(p,'smooth_bands',false,@islogical);
addParameter(p,'export_pdf','',@(s)ischar(s)||isstring(s));
addParameter(p,'tile_shape',[2 3],@(v)isnumeric(v)&&numel(v)==2);
addParameter(p,'date_jitter_days',6,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'y_jitter_frac',0.03,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
parse(p,varargin{:});
args = p.Results;

results_root = args.results_root;
cache_path   = fullfile(results_root,'cache','hbm_cache.mat');
if ~isfile(cache_path)
    error('未找到缓存：%s（请先运行 prep_hbm_inputs）', cache_path);
end
S = load(cache_path); cache = S.cache;

% 目标集合（默认 6 个）
if isempty(args.targets)
    targets = cache.active_targets;
else
    targets = cellstr(string(args.targets));
end
if isempty(targets)
    error('未在缓存中找到 active_targets，或未指定 targets。');
end

% 读取波级信息（features）
Feat = load_features_table(args.features_mat);
if ~istable(Feat)
    error('无法从 %s 读取表变量。', args.features_mat);
end

% 列名标准化与必要字段
Feat.variant = string(strtrim(Feat.variant));
v_levels = ["Original","Delta","Omicron"];
Feat.variant(~ismember(Feat.variant, v_levels)) = "Original";

% 代表时间：优先 start_date / wave_start / t_mid（datetime）
date_col = pick_date_col(Feat);
xdate = Feat.(date_col);
if ~isdatetime(xdate)
    try
        % 常见情形：excel datenum 或字符串
        if isnumeric(xdate)
            xdate = datetime(xdate, 'ConvertFrom','datenum'); %#ok<DATNM>
        else
            xdate = datetime(string(xdate));
        end
    catch
        warning('未找到合法的日期列，将用波索引构造伪日期。');
        xdate = datetime(2020,1,1) + days(double(Feat.wave_idx));
    end
end
Feat.xdate = xdate;

% 点径量纲：扩展候选并级联兜底
metric = lower(string(args.point_size_metric));
mcol = '';
switch metric
    case "cases_total"
        mcol = pick_first_existing(Feat, ["cases_total","total_cases","cases_sum","cases","cum_cases"]);
        if isempty(mcol), mcol = pick_first_existing(Feat, ["duration_days","duration","days"]); end
        if isempty(mcol), mcol = pick_first_existing(Feat, ["w_tail_14","w_tail","w_tail_7","w_tail14","w_tail07"]); end
    case "duration_days"
        mcol = pick_first_existing(Feat, ["duration_days","duration","days"]);
        if isempty(mcol), mcol = pick_first_existing(Feat, ["cases_total","total_cases","cases_sum","cases","cum_cases"]); end
        if isempty(mcol), mcol = pick_first_existing(Feat, ["w_tail_14","w_tail","w_tail_7","w_tail14","w_tail07"]); end
    case "wtail"
        mcol = pick_first_existing(Feat, ["w_tail_14","w_tail","w_tail_7","w_tail14","w_tail07"]);
        if isempty(mcol), mcol = pick_first_existing(Feat, ["cases_total","total_cases","cases_sum","cases","cum_cases"]); end
        if isempty(mcol), mcol = pick_first_existing(Feat, ["duration_days","duration","days"]); end
    otherwise
        mcol = pick_first_existing(Feat, ["cases_total","total_cases","cases_sum","cases","cum_cases"]);
        if isempty(mcol), mcol = pick_first_existing(Feat, ["duration_days","duration","days"]); end
        if isempty(mcol), mcol = pick_first_existing(Feat, ["w_tail_14","w_tail","w_tail_7","w_tail14","w_tail07"]); end
end

if isempty(mcol)
    warning('找不到点径度量列，改用常数大小。');
    Feat.msize = 28 * ones(height(Feat),1);  % 统一中号
else
    xraw = double(Feat.(mcol));
    xraw(~isfinite(xraw) | xraw<0) = NaN;
    % 对数压缩 + 标准化到 [0,1] → marker size（pt）
    xlog = log1p(xraw);
    xr   = normalize_to_01(xlog);
    % 保证至少有区分度
    Feat.msize = 14 + 42* xr;   % 14–56 pt
end

% 颜色（与其它图统一、可改）
cols = struct('Original',[29 120 193]/255, ...   % 蓝
              'Delta',   [233 135 59 ]/255, ...  % 橙
              'Omicron', [111 191 115]/255);     % 绿

% 布局
nT = numel(targets);
nrows = args.tile_shape(1); ncols = args.tile_shape(2);
if nrows*ncols < nT
    nrows = ceil(nT / ncols);
end

% 画布
f = figure('Color','w','Units','pixels','Position',[50 50 1380 820], ...
    'Name','Fig3c Variant timeline');
tl = tiledlayout(f, nrows, ncols, 'Padding','compact','TileSpacing','compact');

annotation(f,'textbox',[0 0 1 0.04], 'String', ...
    "Variant timeline bands + wave scatters", ...
    'FontName','Helvetica','FontWeight','bold','FontSize',16);

% sgtitle(tl, "Variant timeline bands + wave scatters", ...
%     'FontName','Helvetica','FontWeight','bold','FontSize',16);

% 公共 x 轴范围
xmin = min(Feat.xdate); xmax = max(Feat.xdate);
xpad = days(20);
xlim_all_dt = [xmin - xpad, xmax + xpad];
xlim_all_dn = datenum(xlim_all_dt);

% 逐目标作图
for ti = 1:nT
    yname = targets{ti};
    nexttile(tl, ti);
    ax = gca; hold(ax,'on');
    ax.Color = 'w'; ax.FontName = 'Helvetica'; ax.LineWidth = 0.8;
    ax.TickDir = 'out'; grid(ax,'on'); ax.GridAlpha = 0.07; ax.GridColor = [0 0 0];
    ax.YAxis.FontSize = 14; ax.XAxis.FontSize = 14;

    % 读取变异株 RE（常数，不含时间）
    ranef_path = fullfile(results_root, 'M0_main', yname, 'ranef_variant.csv');
    if ~isfile(ranef_path)
        warning('缺少 %s 的 ranef_variant.csv，跳过该目标。', yname);
        title(ax, sprintf('%s  (missing)', yname), 'Interpreter','none');
        continue;
    end
    Tvar = readtable(ranef_path);
    Tvar.variant = string(Tvar.variant);

    % 统一变异株顺序并取统计量
    [~, idxv] = ismember(v_levels, Tvar.variant);
    muV = nan(1,3); loV = nan(1,3); hiV = nan(1,3); medV = nan(1,3);
    for vv = 1:3
        if idxv(vv) > 0
            muV(vv)  = Tvar.mean(idxv(vv));
            loV(vv)  = Tvar.ci2_5(idxv(vv));
            hiV(vv)  = Tvar.ci97_5(idxv(vv));
            medV(vv) = Tvar.median(idxv(vv));
        end
    end

    % 各变异株的时间段（根据波的日期）
    band_spans = NaT(3,2);  % datetime NaT，修复原始 double 分配错误
    for vv = 1:3
        mask = Feat.variant == v_levels(vv);
        if any(mask)
            band_spans(vv,1) = min(Feat.xdate(mask));
            band_spans(vv,2) = max(Feat.xdate(mask));
        end
    end

    % 画“时间带”（CI）+ 中位线（转 datenum 放到 numeric 轴）
    for vv = 1:3
        if ~isfinite(medV(vv)) || any(isnat(band_spans(vv,:)))
            continue;
        end
        x1 = band_spans(vv,1); x2 = band_spans(vv,2);
        if x1==x2, x2 = x1 + days(7); end
        xx_dn = datenum([x1 x2 x2 x1]);
        yy    = [loV(vv) loV(vv) hiV(vv) hiV(vv)];

        % 半透明 CI 带
        patch('XData',xx_dn, 'YData',yy, ...
            'FaceColor',cols.(char(v_levels(vv))), ...
            'FaceAlpha',0.18, 'EdgeColor','none', ...
            'Tag',sprintf('band_CI_%s', v_levels(vv)));

        % 中位线
        plot(datenum([x1 x2]), [medV(vv) medV(vv)],'-', ...
            'Color',cols.(char(v_levels(vv))), 'LineWidth',2.0, ...
            'Tag',sprintf('band_median_%s', v_levels(vv)));
    end

    % y 轴范围
    y_all = [loV hiV medV]; y_all = y_all(isfinite(y_all));
    if isempty(y_all), y_all = 0; end
    ypad = 0.12 * range_or_one(y_all);
    ylim(ax, [min(y_all)-ypad, max(y_all)+ypad]);

    % 叠加波散点
    rng(42 + ti);
    yspan = diff(ylim(ax)); jitter_amp = args.y_jitter_frac * yspan;

    for vv = 1:3
        mask = Feat.variant == v_levels(vv);
        if ~any(mask), continue; end

        xi_dt = Feat.xdate(mask);
        if args.date_jitter_days > 0
            xi_dt = xi_dt + days((rand(size(xi_dt))-0.5)*2*args.date_jitter_days);
        end
        xi_dn = datenum(xi_dt);

        yi0 = medV(vv);
        if ~isfinite(yi0), yi0 = median(y_all,'omitnan'); end
        yi = yi0 + (rand(size(xi_dn))-0.5)*2*jitter_amp;

        si = Feat.msize(mask);

        scatter(xi_dn, yi, si, ...
            'Marker','o','MarkerFaceColor',cols.(char(v_levels(vv))), ...
            'MarkerEdgeColor',[0 0 0], 'LineWidth',0.15, ...
            'MarkerFaceAlpha',0.65, 'Tag',sprintf('wave_%s', v_levels(vv)));
    end

    % 坐标轴样式 & 注记
    xlim(ax, xlim_all_dn);
    datetick(ax,'x','yyyy','keeplimits');  % 可改 'mmm yyyy'
    xlabel(ax, 'Date', 'FontSize',13);
    ylabel(ax, 'Variant random effect (transformed scale)', 'FontSize',13);

    % 零线
    plot(ax, xlim(ax), [0 0],'k-','LineWidth',0.6,'Color',[0 0 0 0.15],'Tag','zero_line');

    % 图例
    draw_compact_legend(ax, cols);

    % 小标题（目标名 + 变换口径）
    tf = cache.tf_info.(yname);
    title(ax, sprintf('%s  (tf: %s)', yname, tf), 'Interpreter','none', ...
        'FontSize',14,'FontWeight','bold');
end

% 空白子图补位
for k = nT+1 : nrows*ncols
    nexttile(tl,k); axis off;
end

% 页脚说明
annotation(f,'textbox',[0 0 1 0.04], 'String', ...
    'Bands: variant-level posterior (median + 95% CI). Dots: waves (size encodes scale). All on transformed scale.', ...
    'EdgeColor','none','HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontName','Helvetica','FontSize',13,'Color',[0 0 0]);

% 导出（可选）
if ~isempty(args.export_pdf)
    try
        set(f,'PaperPositionMode','auto');
        print(f, args.export_pdf,'-dpdf','-painters');  % 矢量
        fprintf('[Fig3c] 已导出 PDF：%s\n', args.export_pdf);
    catch ME
        warning(sprintf('导出 PDF 失败：%s', ME.message));
    end
end
end % ===== 主函数结束 =====


% ==================== 工具函数 ====================

function T = load_features_table(path)
T = [];
if ~isfile(path), return; end
S = load(path);
if isfield(S,'Features') && istable(S.Features)
    T = S.Features; return;
end
fn = fieldnames(S);
for i=1:numel(fn)
    if istable(S.(fn{i})), T = S.(fn{i}); return; end
end
end

function name = pick_date_col(T)
cands = ["start_date","wave_start","t_mid","date","start","begin"];
name = '';
for i=1:numel(cands)
    if ismember(cands(i), T.Properties.VariableNames)
        name = char(cands(i)); return;
    end
end
vn = string(T.Properties.VariableNames);
hit = find(contains(lower(vn),"date") | contains(lower(vn),"start"),1);
if ~isempty(hit), name = char(vn(hit)); else, name = 'wave_idx'; end
end

function nm = pick_first_existing(T, cand)
nm = '';
for i=1:numel(cand)
    if ismember(cand(i), T.Properties.VariableNames)
        nm = char(cand(i)); return;
    end
end
end

function r = normalize_to_01(x)
x = double(x);
if all(~isfinite(x)), r = zeros(size(x)); return; end
mn = min(x,[],'omitnan'); mx = max(x,[],'omitnan');
if ~isfinite(mn) || ~isfinite(mx) || mx==mn
    r = zeros(size(x));
else
    r = (x - mn) ./ (mx - mn);
    r(~isfinite(r)) = 0;
end
end

function w = range_or_one(y)
if isempty(y) || all(~isfinite(y)), w = 1; return; end
w = max(y) - min(y);
if ~(isfinite(w) && w>0), w = 1; end
end

function draw_compact_legend(ax, cols)
xl = xlim(ax); yl = ylim(ax);
dx = xl(2)-xl(1); dy = yl(2)-yl(1);
x0 = xl(1) + 0.985*dx;
y0 = yl(2) - 0.03*dy;

txt = {'Original/Alpha','Delta','Omicron'};
col = {cols.Original, cols.Delta, cols.Omicron};

ystep = 0.08*dy;
for i=1:numel(txt)
    line(ax, [x0-0.05*dx, x0-0.004*dx], [y0 - (i-1)*ystep, y0 - (i-1)*ystep], ...
        'Color', col{i}, 'LineWidth', 5, 'Tag',['legend_' txt{i}]);
    text(ax, x0-0.002*dx, y0 - (i-1)*ystep, txt{i}, ...
        'FontName','Helvetica','FontSize',12, ...
        'HorizontalAlignment','right','VerticalAlignment','middle', ...
        'Color',[0 0 0], 'Tag',['legend_txt_' txt{i}]);
end
end
