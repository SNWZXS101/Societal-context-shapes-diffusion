function uq_smf_features_main()
% UQ/Robustness for SMF features (base-aligned; TeX labels; popup figures)
% 依赖：integrate_seir_odeint_style.m / compute_features_from_wave.m /
%       run_extract_all_features.m 置于路径上；同目录有 waves_dataset.mat

%% ---------------- Config: perturbation grid ----------------
% base 参数与 Methods 一致
BASE = struct( ...
    'name',      "base", ...
    'c_eps',     0.3, 'q_early', 0.30, ...
    'c_eps_hi',  0.3, 'c_eps_lo', 0.3, ...  % 用于构造 ε_hi/ε_neg/ε_nz
    'alpha_hi',  0.2, 'alpha_lo', 0.2, ...
    'eps_mode',  "hybrid-min", ...
    'm_neg', 1,  'm_nz', 1, ...
    'tail_ref', 14);                         % 右尾窗口（天）

% 小工具：复制并批量改字段
makeCfg = @(S, varargin) local_makeCfg(S, varargin{:});

% 扰动集合（与你写法一致）
GRID = [ ...
    BASE, ...
    makeCfg(BASE,'name',"c_0p1",    'c_eps_hi',0.10,'c_eps_lo',0.10), ...
    makeCfg(BASE,'name',"c_0p5",    'c_eps_hi',0.50,'c_eps_lo',0.50), ...
    makeCfg(BASE,'name',"alpha_0p1",'alpha_hi',0.10,'alpha_lo',0.10), ...
    makeCfg(BASE,'name',"alpha_0p3",'alpha_hi',0.30,'alpha_lo',0.30), ...
    makeCfg(BASE,'name',"m0",       'm_neg',0,'m_nz',0), ...
    makeCfg(BASE,'name',"m2",       'm_neg',2,'m_nz',2), ...
    makeCfg(BASE,'name',"m3",       'm_neg',3,'m_nz',3), ...
    makeCfg(BASE,'name',"tail7",    'tail_ref',7), ...
    makeCfg(BASE,'name',"tail10",   'tail_ref',10) ...
];

fprintf('[UQ] Running feature extraction for %d settings...\n', numel(GRID));
Results = struct('name',[],'Features',[]);

for k = 1:numel(GRID)
    A = GRID(k);
    fprintf('  -> %s\n', A.name);
    Ft = run_extract_all_features('waves_dataset.mat', ...
        'outMat','', 'outCsv','', ...
        'useParallel', false, 'verbose', false, ...
        'c_eps',A.c_eps,'q_early',A.q_early, ...
        'c_eps_hi',A.c_eps_hi,'c_eps_lo',A.c_eps_lo, ...
        'alpha_hi',A.alpha_hi,'alpha_lo',A.alpha_lo, ...
        'eps_mode',A.eps_mode, ...
        'm_neg',A.m_neg,'m_nz',A.m_nz, ...
        'tail_ref',A.tail_ref);                 % ← 右尾窗口下传

    if ~ismember('wave_id', Ft.Properties.VariableNames)
        error('Features 缺少 wave_id，无法对齐比较。');
    end
    Ft.wave_id = string(Ft.wave_id);
    Results(k).name     = A.name;
    Results(k).Features = Ft;
end

% 找 base（string 安全）
names = arrayfun(@(r) string(r.name), Results);
idx_base = find(names=="base", 1);
if isempty(idx_base)
    error('未找到名为 "base" 的基线设置。可用: %s', strjoin(names, ', '));
end
Runs = names; M = numel(Runs);

% 以 base 顺序对齐
BaseT = Results(idx_base).Features;
WIDs  = BaseT.wave_id;
for k = 1:numel(Results)
    T = Results(k).Features; [tf, loc] = ismember(WIDs, T.wave_id);
    T2 = table(); T2.wave_id = WIDs;
    for j = 1:numel(T.Properties.VariableNames)
        vj = string(T.Properties.VariableNames{j});
        if vj ~= "wave_id"
            col = nan(height(T2),1);
            col(tf) = toCol(T(loc(tf), vj));
            T2.(vj) = col;
        end
    end
    Results(k).Features = T2;
end

% ---------- 派生要用的变量清单（存在才取） ----------
allCols = string(BaseT.Properties.VariableNames);
wanted = string({ ...
    't_cross','kappa_cross','delta_early', ...
    'A_plus','A_minus','A_minus_eps','rho_eps', ...
    'tf1_mean_early','tf2_peak', ...
    'w_tail_14', ...
    'flag_multipart'});
KEYVARS   = wanted(ismember(wanted, allCols));
wanted_d  = string({'t_cross','kappa_cross','delta_early','A_plus','rho_eps'});
DIST_VARS = wanted_d(ismember(wanted_d, KEYVARS));
if isempty(KEYVARS)
    isNumCol = varfun(@(x) isnumeric(x) || islogical(x), BaseT, 'OutputFormat','uniform');
    KEYVARS = string(BaseT.Properties.VariableNames(isNumCol));
end
if isempty(DIST_VARS), DIST_VARS = KEYVARS(1:min(5,numel(KEYVARS))); end

% ---------- 计算稳健性量 ----------
K = numel(KEYVARS);
Rho = nan(K, M); Dmed = nan(K, M); PctMiss = nan(K, M);
for m = 1:M
    Tm = Results(m).Features;
    for i = 1:K
        v = KEYVARS(i);
        xb = pickVar(BaseT, v); xm = pickVar(Tm, v);
        mask = isfinite(xb) & isfinite(xm);
        if any(mask)
            Rho(i,m) = corr(xb(mask), xm(mask), 'Type','Spearman','Rows','pairwise');
            Dmed(i,m)= median(xm(mask)) - median(xb(mask));
        end
        PctMiss(i,m) = mean(~isfinite(xm));
    end
end

% ---------- FIG S2: 分布稳定性 ----------
figure('Name','S2-Distributions','Color','w');
tl = tiledlayout(numel(DIST_VARS), 1, 'TileSpacing','compact','Padding','compact');
prettyRuns = arrayfun(@pretty_run_label, Runs, 'UniformOutput', false);
for r = 1:numel(DIST_VARS)
    nexttile; hold on;
    v = DIST_VARS(r);
    X = []; G = [];
    for m = 1:M
        xm = pickVar(Results(m).Features, v); xm = xm(isfinite(xm));
        X = [X; xm(:)]; G = [G; repmat(m, numel(xm), 1)];
    end
    if isempty(X)
        title(sprintf('%s (no data)', v), 'Interpreter','tex','FontSize',14);
        axis off; continue;
    end
    boxchart(G, X, 'BoxFaceAlpha',0.35, 'BoxWidth',0.6);
    jitter = 0.25*(rand(size(X))-0.5);
    scatter(double(G)+jitter, X, 8, 'filled', ...
        'MarkerFaceAlpha',0.18, 'MarkerEdgeAlpha',0.12);
    xlim([0.5 M+0.5]);
    set(gca,'XTick',1:M,'XTickLabel',prettyRuns,'XTickLabelRotation',30, ...
        'TickLabelInterpreter','tex','FontSize',14);
    ylabel(pretty_var_label(v), 'Interpreter','tex','FontSize',16); grid on;
end
title(tl, 'Key features: distribution stability across perturbations', ...
    'FontSize',16,'Interpreter','tex');

% ---------- FIG S3: Spearman 热图（带数值与缺失星标） ----------
figure('Name','S3-SpearmanHeat','Color','w');
K2 = numel(KEYVARS);
imagesc(Rho, [0 1]); axis tight; colormap(parula);
cb = colorbar; cb.Label.String = '\rho_{Spearman}'; cb.Label.Interpreter = 'tex';
set(gca,'YTick',1:K2,'YTickLabel',arrayfun(@pretty_var_label, KEYVARS, 'UniformOutput', false), ...
        'XTick',1:M,'XTickLabel',prettyRuns, ...
        'TickLabelInterpreter','tex','FontSize',14,'XTickLabelRotation',30);
title('Spearman rank correlation vs base','Interpreter','tex','FontSize',16);
grid on; set(gca,'Layer','top'); hold on;
for i = 1:K2
    for j = 1:M
        val = Rho(i,j); miss = PctMiss(i,j);
        if isfinite(val)
            txt = sprintf('%.2f', val); col = [0 0 0]; if val >= 0.60, col = [1 1 1]; end
        else
            txt = '–'; col = [0.3 0.3 0.3];
        end
        text(j, i, txt, 'HorizontalAlignment','center', ...
            'FontSize',12, 'FontWeight','bold', 'Color', col, 'Interpreter','tex');
        if isfinite(miss) && miss > 0.10
            text(j+0.26, i+0.22, '*', 'HorizontalAlignment','center', ...
                'FontSize',16, 'FontWeight','bold', 'Color','k');
        end
    end
end
hold off

% ---------- FIG S4: Crossing 组成 + 多波占比 ----------
figure('Name','S4-Outcomes+Multipart','Color','w');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% (a) crossing outcome 组成
nexttile; 
Tprop = zeros(3, M); % [neg; nz; censored]
for m = 1:M
    em = pickVar(Results(m).Features, "event_observed");
    tm = pickVar(Results(m).Features, "event_type");
    Tprop(1,m) = mean(isfinite(em) & em>0.5 & isfinite(tm) & (tm>0.5));      % neg
    Tprop(2,m) = mean(isfinite(em) & em>0.5 & isfinite(tm) & (abs(tm)<0.5)); % nz
    Tprop(3,m) = mean(isfinite(em) & (em<0.5));                               % censored
end
bar((1:M), Tprop', 'stacked'); ylim([0 1]); grid on;
set(gca,'XTick',1:M,'XTickLabel',prettyRuns,'XTickLabelRotation',30, ...
    'TickLabelInterpreter','tex','FontSize',14);
legend({'Negative phase','Near-zero','Censored'}, 'Location','eastoutside','FontSize',12);
ylabel('Proportion','FontSize',16); 
title('Crossing outcome composition across perturbations','FontSize',16,'Interpreter','tex');

% (b) 多波占比
nexttile;
rate = nan(1,M);
for m = 1:M
    fm = pickVar(Results(m).Features, "flag_multipart");
    rate(m) = mean(isfinite(fm) & fm>0.5);
end
bar(1:M, rate); ylim([0 1]); grid on;
set(gca,'XTick',1:M,'XTickLabel',prettyRuns,'XTickLabelRotation',30, ...
    'TickLabelInterpreter','tex','FontSize',14);
ylabel('Share of multipart waves','FontSize',16,'Interpreter','tex');
title('Multipart-flag rate across perturbations','FontSize',16,'Interpreter','tex');

fprintf('[UQ] Done. Figures popped. No files saved.\n');
end

%% ---------------- helpers ----------------
function S = local_makeCfg(S, varargin)
n = numel(varargin);
if mod(n,2)~=0, error('makeCfg 需要成对的 field/value 参数'); end
for i = 1:2:n
    fld = varargin{i}; val = varargin{i+1};
    S.(char(fld)) = val;
end
end

function v = toCol(x)
try
    if istable(x)
        if width(x)==1, v = double(table2array(x)); else, v = double(x{:,1}); end
    else, v = double(x);
    end
catch
    try, v = str2double(string(x)); catch, v = nan(size(x,1),1); end
end
v = v(:);
end

function x = pickVar(T, varname)
vn = string(T.Properties.VariableNames);
idx = find(vn==string(varname),1);
if isempty(idx), x = nan(height(T),1); else, x = toCol(T(:,idx)); end
end

function lbl = pretty_var_label(v)
v = string(v);
switch v
    case "t_cross",           lbl = 't_{cross}';
    case "kappa_cross",       lbl = '\kappa_{cross}';
    case "delta_early",       lbl = '\Delta_{early}';
    case "delta_early_rel",   lbl = '\Delta_{early}/T';
    case "A_plus",            lbl = 'A^{+}';
    case "A_minus",           lbl = 'A^{-}';
    case "A_minus_eps",       lbl = 'A^{-}_{\epsilon}';
    case "rho_eps",           lbl = '\rho_{\epsilon}';
    case "tf1_mean_early",    lbl = 'f_{1}^{early}';
    case "tf2_peak",          lbl = 'f_{2}^{peak}';
    case "w_tail_14",         lbl = 'w_{tail,14}';
    case "w_tail_rel",        lbl = 'w_{tail,rel}';
    case "flag_multipart",    lbl = 'flag_{multipart}';
    case "event_time",        lbl = 't_{event}';
    case "event_observed",    lbl = '\mathrm{obs}';
    case "event_type",        lbl = '\mathrm{type}';
    otherwise,                lbl = v;
end
end

function lbl = pretty_run_label(name)
s = string(name);
switch s
    case "base",       lbl = 'base';
    case "c_0p1",      lbl = 'c_{hi}=c_{lo}=0.10';
    case "c_0p5",      lbl = 'c_{hi}=c_{lo}=0.50';
    case "alpha_0p1",  lbl = '\alpha_{hi}=\alpha_{lo}=0.10';
    case "alpha_0p3",  lbl = '\alpha_{hi}=\alpha_{lo}=0.30';
    case "m0",         lbl = 'm_{nz}=m_{neg}=0';
    case "m2",         lbl = 'm_{nz}=m_{neg}=2';
    case "m3",         lbl = 'm_{nz}=m_{neg}=3';
    case "tail7",      lbl = 'w_{tail}=7';
    case "tail10",     lbl = 'w_{tail}=10';
    otherwise,         lbl = s;
end
end
