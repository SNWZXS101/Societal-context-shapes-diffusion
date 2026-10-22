function Plot_PPC_and_Diagnostics(root, topK)
% 左：PPC（变换后尺度）；右：R-hat/ESS 诊断（从 summary.csv 的 rhat/ess 读取）
% 用法：Plot_PPC_and_Diagnostics('results_hbm/M0_main', 12)

if nargin<1, root = 'results_hbm/M0_main'; end
if nargin<2, topK = 12; end

out_dir = fullfile(root, '..', 'figs_ppc_diag');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
nature_style();

% 读 tf_info 作为备选
cache_path = fullfile(fileparts(root), 'cache', 'hbm_cache.mat');
tf_info = struct();
if exist(cache_path,'file')
    S = load(cache_path); 
    if isfield(S,'cache') && isfield(S.cache,'tf_info')
        tf_info = S.cache.tf_info;
    end
end

ds = dir(root); ds = ds([ds.isdir]);
targets = setdiff({ds.name},{'.','..'});

for i=1:numel(targets)
    yname = targets{i};
    ydir  = fullfile(root, yname);
    f_ppc = fullfile(ydir,'ppc.csv');
    f_sum = fullfile(ydir,'summary.csv');
    if ~exist(f_ppc,'file') || ~exist(f_sum,'file')
        fprintf('[skip] %s 缺 ppc/summary\n', yname); continue;
    end

    PPC = readtable(f_ppc);
    SUM = readtable(f_sum);

    % 取 tf
    tf = '';
    if ismember('tf', PPC.Properties.VariableNames)
        tf = char(string(PPC.tf(1)));
    end
    if isempty(tf)
        fn = fieldnames(tf_info);
        hit = find(strcmpi(fn, yname),1);
        if ~isempty(hit), tf = tf_info.(fn{hit}); else, tf = 'identity'; end
    end
    tflabel = make_tf_label(yname, tf);

    % ---------- 画布 ----------
    fig = figure('Color','w','Position',[100 100 1100 520]);
    TL  = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

    % ===== 左：PPC =====
    axP = nexttile(TL,1); hold(axP,'on');
    x  = PPC.y_obs;  y = PPC.yhat_p50;
    lo = y - PPC.yhat_p10; hi = PPC.yhat_p90 - y;
    errorbar(axP, x, y, lo, hi, 'LineStyle','none','CapSize',0,'Color',[0 0 0]);
    plot(axP, x, y, 'o', 'MarkerFaceColor',[0 0.45 0.74], 'MarkerEdgeColor','none', 'MarkerSize',5);

    mm = nanminmax([x; y; PPC.yhat_p10; PPC.yhat_p90]);
    plot(axP, mm, mm, 'k--');
    xlim(axP, mm); ylim(axP, mm);

    cover80 = mean(PPC.y_obs >= PPC.yhat_p10 & PPC.y_obs <= PPC.yhat_p90, 'omitnan');
    rmse    = sqrt(mean((y - x).^2, 'omitnan'));
    title(axP, sprintf('Posterior predictive check — %s', tflabel),'Interpreter','none');
    xlabel(axP, sprintf('Observed — %s', tflabel), 'Interpreter','none');
    ylabel(axP, sprintf('Predicted (median) — %s', tflabel), 'Interpreter','none');
    txt = sprintf('80%% coverage = %.2f | RMSE = %.3f | N = %d', cover80, rmse, height(PPC));
    text(axP, 0.02,0.98, txt, 'Units','normalized','VerticalAlignment','top');

    % ===== 右：R-hat / ESS =====
    T2 = tiledlayout(TL,2,1,'TileSpacing','compact','Padding','compact'); T2.Layout.Tile = 2;

    labs_all  = to_cellstr(SUM.param);
    base_names= {'alpha','sigma','sigma_prov','sigma_var','nu'};
    base_rows = ismember(lower(labs_all), lower(base_names));   % 大小写无关
    bet_rows  = startsWith_cell_icase(labs_all, 'beta:');

    SUMbeta   = SUM(bet_rows,:);
    medabs    = abs(SUMbeta.median); medabs(~isfinite(medabs)) = -inf;
    [~,ord]   = sort(medabs,'descend');
    pickK     = min(topK, height(SUMbeta));
    SUMpick   = [SUM(base_rows,:); SUMbeta(ord(1:pickK),:)];

    % 如果没有 rhat/ess 列，补 NaN 并提示
    if ~ismember('rhat', SUMpick.Properties.VariableNames), SUMpick.rhat(:) = NaN; end
    if ~ismember('ess',  SUMpick.Properties.VariableNames), SUMpick.ess(:)  = NaN; end

    labs = to_cellstr(SUMpick.param); yidx = 1:numel(labs);

    % --- R-hat ---
    ax1 = nexttile(T2,1); hold(ax1,'on');
    rh = SUMpick.rhat; plot(ax1, rh, yidx, 'o', 'MarkerFaceColor',[0 0.45 0.74],'MarkerEdgeColor','none','MarkerSize',5);
    xline(ax1, 1.01, 'r--', '0.01 threshold');
    set(ax1,'YTick',yidx,'YTickLabel',labs,'YDir','reverse'); grid(ax1,'on'); ax1.GridAlpha=0.2;
    if all(~isfinite(rh)), xlim(ax1,[0.98 1.04]); text(ax1,0.02,0.1,'(No rhat in summary.csv)','Units','normalized'); end
    title(ax1,'Split R-hat (↓ 好)','FontWeight','normal');

    % --- ESS ---
    ax2 = nexttile(T2,2); hold(ax2,'on');
    ess = SUMpick.ess; plot(ax2, ess, yidx, 's', 'MarkerFaceColor',[0 0.45 0.74],'MarkerEdgeColor','none','MarkerSize',5);
    xline(ax2, 400, 'r--', 'ESS=400');
    set(ax2,'YTick',yidx,'YTickLabel',labs,'YDir','reverse'); grid(ax2,'on'); ax2.GridAlpha=0.2;
    if all(~isfinite(ess)), xlim(ax2,[0 1000]); text(ax2,0.02,0.1,'(No ESS in summary.csv)','Units','normalized'); end
    xlabel(ax2,'ESS (↑ 好)'); title(ax2,'Effective sample size','FontWeight','normal');

    % 保存
    png_path = fullfile(out_dir, sprintf('%s_ppc_diagnostics.png', yname));
    pdf_path = fullfile(out_dir, sprintf('%s_ppc_diagnostics.pdf', yname));
    try
        exportgraphics(fig, png_path, 'Resolution', 200);
        exportgraphics(fig, pdf_path, 'ContentType','vector');
    catch
        set(fig,'PaperPositionMode','auto');
        print(fig, png_path, '-dpng','-r200');
        print(fig, pdf_path, '-dpdf');
    end
%     close(fig);
end

fprintf('[PPC+Diag] 输出到：%s\n', out_dir);

% ---------- 工具 ----------
function lbl = make_tf_label(yname, tf)
nm = char(string(yname));
switch lower(tf)
    case 'log',    lbl = sprintf('log(%s)', nm);
    case 'abslog', lbl = sprintf('log|%s|', nm);
    case 'asinh',  lbl = sprintf('asinh(%s)', nm);
    otherwise,     lbl = nm;
end
end

function c = to_cellstr(x)
if iscellstr(x), c = x;
elseif isstring(x), c = cellstr(x);
elseif iscell(x)
    c = cell(size(x));
    for ii=1:numel(x), c{ii} = char(x{ii}); end
else
    c = cellstr(x);
end
end

function tf = startsWith_cell_icase(ca, prefix)
tf = false(size(ca));
for ii=1:numel(ca)
    s = lower(ca{ii}); tf(ii) = strncmp(s, lower(prefix), length(prefix));
end
end

function mm = nanminmax(v)
v = v(:); v = v(isfinite(v));
if isempty(v), mm = [0 1]; else, mm = [floor(min(v)), ceil(max(v))]; end
end

function nature_style()
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
set(0,'DefaultAxesFontSize',10); set(0,'DefaultLineLineWidth',1.25);
set(0,'DefaultAxesBox','off'); set(0,'DefaultFigurePaperPositionMode','auto');
end
end
