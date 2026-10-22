function H = plot_mechanistic_ridge_minmax1(varargin)
% Mechanistic Ridge (min–max quintiles, per-panel autoscale)
% 新增：导出 (Qmax vs Qmin) 的时间位移与幅度比到 CSV
%
% Key params:
%   'features_aug'  : 'features_aug.mat' (table Features)
%   'features_raw'  : 'features.mat'     (cell f_cell, optional t_cell)
%   'covariates'    : cellstr of covariate names (2–6 推荐)
%   'variants'      : {'Alpha','Delta','Omicron'}
%   'nGrid'         : 200
%   'ampNorm'       : 'Aplus' | 'fmax'   (绘图使用的幅度归一)
%   'legendOn'      : true
%   'renderer'      : 'painters'
%   'out_csv'       : 'ridge_contrasts.csv'   % <— 新增：导出文件名
%
% 输出: H.fig, H.ax (R×C), H.meta, H.contrasts_table

%% ---- args
p = inputParser;
p.addParameter('features_aug','features_aug.mat');
p.addParameter('features_raw','features.mat');
p.addParameter('covariates',{'urban rate','public transit volume','disposable income_pc'});
p.addParameter('variants',{'Alpha','Delta','Omicron'});
p.addParameter('nGrid',200);
p.addParameter('ampNorm','Aplus');      % 绘图用
p.addParameter('legendOn',true);
p.addParameter('renderer','painters');
p.addParameter('out_csv','ridge_contrasts.csv');  % 新增
p.parse(varargin{:});
A = p.Results;
covars   = cellstr(string(A.covariates));
variants = canon_variant(cellstr(string(A.variants)));
ng       = A.nGrid;

%% ---- colors
try
    cold = smf_colorlib('mainA8'); cold = cold(7,:);   % Qmin
    warm = smf_colorlib('mainA8'); warm = warm(1,:);   % Qmax
catch
    L = lines(8); cold = L(2,:); warm = L(6,:);
end
Cline = [0.1 0.1 0.1];

%% ---- load Features & f_cell
T = load_features_table(A.features_aug);
[colProv,colVar] = pick_cols(T, ["province","Province"], ["variant","Variant"]);
Tvariant = canon_variant(cellstr(string(T.(colVar))));
colWave  = pick_one(T,["wave_idx","wave","waveId","waveid"]);

if strlength(colWave)==0
    idxWave = (1:height(T))';
else
    idxWave = double(T.(colWave));
end

Sraw   = load(A.features_raw);
f_cell = pick_field(Sraw, {'f_cell','f','fcurve','f_curves','fcell','f_vals'});
t_cell = [];
if isfield(Sraw,'t_cell'), t_cell = Sraw.t_cell; end
assert(iscell(f_cell),'f_cell must be a cell array');

get_curve = @(i) take_f(f_cell, t_cell, idxWave(i), i);  % returns [t,f] on native grid

%% ---- build unified [0,1] grids (两套：按绘图开关、以及固定 Aplus 版用于统计)
tg   = linspace(0,1,ng);
Yg   = nan(height(T), ng);       % 按 A.ampNorm
YgAp = nan(height(T), ng);       % 固定 Aplus（用于对比指标）

for i=1:height(T)
    [~,f0] = get_curve(i);
    if isempty(f0), continue; end
    t0 = linspace(0,1,numel(f0));

    % --- A.ampNorm（绘图）
    if strcmpi(A.ampNorm,'fmax')
        denom = max(abs(f0));
    else
        denom = trapz(t0, max(f0,0));
    end
    if ~isfinite(denom) || denom<=eps, denom = 1; end
    fN  = f0 / denom;
    Yg(i,:) = interp1(t0, fN, tg, 'makima', 'extrap');

    % --- Aplus（用于统计对比，一致口径）
    denomA = trapz(t0, max(f0,0));
    if ~isfinite(denomA) || denomA<=eps, denomA = 1; end
    fA  = f0 / denomA;
    YgAp(i,:) = interp1(t0, fA, tg, 'makima', 'extrap');
end

%% ---- layout (绘图与之前一致)
R = numel(covars); C = numel(variants);
fig = figure('Color','w','Name','Mechanistic Ridge (min–max)','Renderer',char(A.renderer));
tlo = tiledlayout(fig, R, C, 'TileSpacing','compact','Padding','compact');
ax  = gobjects(R,C);

for j=1:C
    pos = tilebox(nexttile(tlo,(1-1)*C+j), fig, [0 0 0 0]); % 第一行第 j 列位置
    annotation(fig,'textbox',[pos(1) pos(2)+pos(4)+0.01 pos(3) 0.035], ...
        'String', variants{j}, 'Interpreter','none', ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontName','Arial','FontSize',14,'FontWeight','normal','EdgeColor','none');
end

%% ---- pre-allocate contrast collectors
rows_cov  = cell(R*C,1);
rows_var  = cell(R*C,1);
nQmin_arr = zeros(R*C,1);
nQmax_arr = zeros(R*C,1);
t0_min_arr= nan(R*C,1); t0_max_arr = nan(R*C,1); dt_zero_arr = nan(R*C,1);
tc_min_arr= nan(R*C,1); tc_max_arr = nan(R*C,1); dt_cent_arr = nan(R*C,1);
fmax_minA = nan(R*C,1); fmax_maxA = nan(R*C,1); ratioAmp  = nan(R*C,1);
rowPtr = 0;

%% ---- draw each cell (autoscale Y per panel) + 计算对比分数
Lh = struct('min',gobjects(0), 'max',gobjects(0));
for r=1:R
    xcov = pull_covar(T, covars{r});
    ylab = wraplabel(pretty_name(covars{r}), 14);
    for j=1:C
        ax(r,j) = nexttile(tlo, (r-1)*C + j); hold(ax(r,j),'on'); box(ax(r,j),'off');
        set(ax(r,j),'FontName','Arial','TickDir','out','LineWidth',0.8,'Layer','top');
        if j==1
            ylabel(ax(r,j), ylab, 'Interpreter','tex');
        else
            ax(r,j).YColor = 'none';  % 隐藏 y 轴
            ax(r,j).YTick  = [];
        end
        if r==R
            xlabel(ax(r,j),'Normalized time');
        else
            ax(r,j).XColor = 'none';  % 隐藏 x 轴
            ax(r,j).XTick  = [];
        end
        plot(ax(r,j), [0 1], [0 0], '-', 'Color',[0.87 0.87 0.87], 'LineWidth',0.8,'HitTest','off');

        % —— 本列的变异株掩码
        mV = strcmpi(Tvariant, variants{j});
        xv = xcov; xv(~mV) = NaN;
        if all(~isfinite(xv))
            text(ax(r,j),0.5,0.5,'no data','Units','normalized'); 
            continue; 
        end

        % —— quintile edges
        q = quantile_strict(xv, 5);

        % —— pick min & max bins (索引，同步用于两套网格)
        [idxMin,idxMax,nMin,nMax, qminTxt,qmaxTxt] = pick_minmax(xv, q);

        % —— 绘图（按 A.ampNorm）
        ym = rows_mean(Yg(idxMin,:)); 
        yM = rows_mean(Yg(idxMax,:));
        ymax = max(abs([ym yM]), [], 'all', 'omitnan');
        if ~isfinite(ymax) || ymax<=0, ymax = 1; end
        ylim(ax(r,j), 1.05*[-ymax ymax]);
        xlim(ax(r,j), [0 1]);

        [p1,p2,l1] = draw_signed(ax(r,j), tg, ym, cold, Cline);
        [p3,p4,l2] = draw_signed(ax(r,j), tg, yM, warm, Cline);
        uistack([p1 p2 p3 p4],'bottom');
        Lh.min(end+1,1) = l1; %#ok<AGROW>
        Lh.max(end+1,1) = l2;

        % 角落标签
        if nMin>0
            text(ax(r,j), 0.02, 0.98, sprintf('Q_{min} %s (n=%d)',qminTxt,nMin), ...
                'Units','normalized','Color',cold*0.85,'VerticalAlignment','top');
        end
        if nMax>0
            text(ax(r,j), 0.02, 0.02, sprintf('Q_{max} %s (n=%d)',qmaxTxt,nMax), ...
                'Units','normalized','Color',warm*0.85,'VerticalAlignment','bottom');
        end

        % ========= 统计对比（固定 Aplus 口径） =========
        ymA = rows_mean(YgAp(idxMin,:)); 
        yMA = rows_mean(YgAp(idxMax,:));

        % 1) 主峰后首次过零
        t0_min = zero_cross_after_peak(tg, ymA);
        t0_max = zero_cross_after_peak(tg, yMA);
        dt0    = t0_max - t0_min;

        % 2) 正瓣质心（∫t·max(f,0) / ∫max(f,0)）
        tc_min = pos_centroid(tg, ymA);
        tc_max = pos_centroid(tg, yMA);
        dtc    = tc_max - tc_min;

        % 3) 峰值比（Aplus 归一后）
        fmax_qmin = max(ymA);
        fmax_qmax = max(yMA);
        if ~isfinite(fmax_qmin) || fmax_qmin==0
            rat = NaN;
        else
            rat = fmax_qmax / fmax_qmin;
        end

        % —— 收集
        rowPtr = rowPtr + 1;
        rows_cov{rowPtr,1}  = covars{r};
        rows_var{rowPtr,1}  = variants{j};
        nQmin_arr(rowPtr,1) = nMin;
        nQmax_arr(rowPtr,1) = nMax;
        t0_min_arr(rowPtr,1)= t0_min;
        t0_max_arr(rowPtr,1)= t0_max;
        dt_zero_arr(rowPtr,1)= dt0;
        tc_min_arr(rowPtr,1)= tc_min;
        tc_max_arr(rowPtr,1)= tc_max;
        dt_cent_arr(rowPtr,1)= dtc;
        fmax_minA(rowPtr,1) = fmax_qmin;
        fmax_maxA(rowPtr,1) = fmax_qmax;
        ratioAmp(rowPtr,1)  = rat;
    end
end

% 图例
if A.legendOn && ~isempty(Lh.min) && ~isempty(Lh.max)
    lg = legend(ax(end,end), [Lh.min(1) Lh.max(1)], ...
        {'Qmin (lowest quintile available)','Qmax (highest quintile available)'}, ...
        'Orientation','horizontal','Location','southoutside','Box','off');
    lg.ItemTokenSize = [18 9];
end

%% ---- 组装表格、打印并导出 CSV
valid = (1:rowPtr).';
Tcontrast = table( ...
    string(rows_cov(valid)), ...
    string(rows_var(valid)), ...
    nQmin_arr(valid), nQmax_arr(valid), ...
    t0_min_arr(valid), t0_max_arr(valid), dt_zero_arr(valid), ...
    tc_min_arr(valid), tc_max_arr(valid), dt_cent_arr(valid), ...
    fmax_minA(valid), fmax_maxA(valid), ratioAmp(valid), ...
    'VariableNames', {'covariate','variant','n_Qmin','n_Qmax', ...
                      't_zero_Qmin','t_zero_Qmax','dt_zero', ...
                      't_centpos_Qmin','t_centpos_Qmax','dt_centpos', ...
                      'fmax_Aplus_Qmin','fmax_Aplus_Qmax','amp_ratio_Aplus'});

% 打印
fprintf('\n[Mechanistic ridge — min–max contrasts (A+ normalization for metrics)]\n');
for ii = 1:height(Tcontrast)
    fprintf('%-22s | %-8s | n=(%2d,%2d) | dt_zero=% .3f | dt_cent=% .3f | amp_ratio(A+)=% .3f\n', ...
        Tcontrast.covariate(ii), Tcontrast.variant(ii), ...
        Tcontrast.n_Qmin(ii), Tcontrast.n_Qmax(ii), ...
        Tcontrast.dt_zero(ii), Tcontrast.dt_centpos(ii), ...
        Tcontrast.amp_ratio_Aplus(ii));
end
try
    writetable(Tcontrast, A.out_csv);
    fprintf('=> Saved contrasts to %s\n', A.out_csv);
catch ME
    warning(sprintf('Failed to write CSV: %s', ME.message));
end

H = struct('fig',fig,'ax',ax, ...
    'meta',struct('variants',{variants},'covariates',{covars},'ampNorm',A.ampNorm), ...
    'contrasts_table',Tcontrast);

end

%% ================= helpers =================
function T = load_features_table(fpath)
S = load(fpath);
T = [];
fn = fieldnames(S);
k = 0;
while isempty(T) && k < numel(fn)
    k = k + 1;
    v = S.(fn{k});
    if istable(v), T = v; end
end
if isempty(T), error('No table found in %s', fpath); end
end

function [cp,cv] = pick_cols(T, candP, candV)
cp = pick_one(T,candP); cv = pick_one(T,candV);
if strlength(cp)==0||strlength(cv)==0, error('Missing province/variant columns'); end
end

function nm = pick_one(T,cands)
vn = string(T.Properties.VariableNames); nm="";
i = 0;
while i < numel(cands) && strlength(nm)==0
    i = i + 1;
    j = find(strcmpi(vn,cands(i)),1);
    if ~isempty(j), nm = vn(j); end
end
end

function v = canon_variant(v)
v = string(v);
i = 0;
while i < numel(v)
    i = i + 1;
    s = lower(strtrim(v(i)));
    if strcmp(s,'alpha') || any(strcmp(s,["original","ancestral","wuhan","wildtype","wild-type","wild type"]))
        v(i) = "Alpha";
    elseif strcmp(s,'delta')
        v(i) = "Delta";
    elseif strcmp(s,'omicron')
        v(i) = "Omicron";
    else
        v(i) = upper(extractBefore(s,1))+extractAfter(s,1);
    end
end
end

function val = pick_field(S, names)
val = [];
i = 0;
while isempty(val) && i < numel(names)
    i = i + 1;
    if isfield(S, names{i}), val = S.(names{i}); end
end
end

function [t,f] = take_f(fcell, tcell, waveId, rowId)
t=[]; f=[]; idx=[];
if isfinite(waveId) && waveId>=1 && waveId<=numel(fcell), idx = waveId;
elseif rowId>=1 && rowId<=numel(fcell), idx = rowId;
else, return;
end
item = fcell{idx};
if isstruct(item)
    if isfield(item,'f')
        f = double(item.f(:)').';
        f = f(:)';  % row
    elseif isfield(item,'y')
        f = double(item.y(:)').';
        f = f(:)'; 
    else
        f = double(item(:)'); 
    end
else
    f = double(item(:)'); 
end
if isempty(tcell)
    t = linspace(0,1,numel(f));
else
    if numel(tcell)>=idx && ~isempty(tcell{idx})
        t = double(tcell{idx}(:))';
    else
        t = linspace(0,1,numel(f));
    end
end
end

function x = pull_covar(T, nameLike)
vn = string(T.Properties.VariableNames);
key = lower(regexprep(string(nameLike),'[^a-z0-9]+',''));
x = nan(height(T),1);
i = 0;
while i < numel(vn)
    i = i + 1;
    kk = lower(regexprep(vn(i),'[^a-z0-9]+',''));
    if strcmp(kk,key), x = double(T.(vn(i))); return; end
end
i = 0;
while i < numel(vn)
    i = i + 1;
    kk = lower(regexprep(vn(i),'[^a-z0-9]+',''));
    if contains(kk,key), x = double(T.(vn(i))); return; end
end
warning('Covariate not found: %s', string(nameLike));
end

function q = quantile_strict(x, Q)
x = x(isfinite(x));
if numel(x)<2 || all(x==x(1))
    lo = min(x); hi = max(x);
    if ~isfinite(lo) || ~isfinite(hi), lo=0; hi=1; end
    q = linspace(lo, hi + max(1e-9,1e-12*max(1,abs(hi))), Q+1);
    return
end
q = quantile(x, linspace(0,1,Q+1)); q = q(:).';
k = 2;
while k <= numel(q)
    if q(k) <= q(k-1)
        epsJ = max(1e-9, 1e-12*max(1,abs(q(k-1))));
        q(k) = q(k-1) + epsJ;
    end
    k = k + 1;
end
end

function [idxMin,idxMax,nMin,nMax, qminTxt,qmaxTxt] = pick_minmax(x, edges)
bins = discretize(x, edges, 'IncludedEdge','right');
avail = unique(bins(isfinite(bins)));
if isempty(avail)
    idxMin = false(size(x)); idxMax = false(size(x));
    nMin=0; nMax=0; qminTxt=''; qmaxTxt=''; return;
end
qmin = min(avail); qmax = max(avail);
idxMin = (bins==qmin);  nMin = nnz(idxMin);
idxMax = (bins==qmax);  nMax = nnz(idxMax);
qminTxt = sprintf('(Q%d)', qmin);
qmaxTxt = sprintf('(Q%d)', qmax);
end

function y = rows_mean(M)
if isempty(M)
    y = nan(1,size(M,2));
else
    y = nanmean(M,1);
end
end

function [p1,p2,l] = draw_signed(ax, x, y, col, colLine)
if all(~isfinite(y)), p1=gobjects(1); p2=gobjects(1); l=gobjects(1); return; end
y  = double(y);
yp = max(y,0); yn = min(y,0);
p1 = patch(ax, [x fliplr(x)], [zeros(size(x)) fliplr(yp)], col, ...
    'FaceAlpha',0.50,'EdgeColor','none','PickableParts','all','HitTest','on');
p2 = patch(ax, [x fliplr(x)], [zeros(size(x)) fliplr(yn)], col, ...
    'FaceAlpha',0.50,'EdgeColor','none','PickableParts','all','HitTest','on');
l  = plot(ax, x, y, '-', 'Color',colLine, 'LineWidth',1.5, 'PickableParts','all','HitTest','on');
end

function s = pretty_name(nm)
s = regexprep(string(nm),'_',' ');
s = regexprep(s,'(\w)(\w*)','${upper($1)}${$2}');
end

function s = wraplabel(txt, maxch)
if nargin < 2, maxch = 14; end
txt = string(txt);
txt = replace(txt,"_","\_"); txt = replace(txt,"{","\{"); txt = replace(txt,"}","\}");
toks = split(txt); toks = toks(strlength(toks) > 0);
lines = strings(0,1); cur = "";
i = 0;
while i < numel(toks)
    i = i + 1; t = toks(i);
    if strlength(cur) == 0, cand = t; else, cand = cur + " " + t; end
    if strlength(cand) <= maxch, cur = cand;
    else
        if strlength(cur) > 0, lines(end+1,1) = cur; end %#ok<AGROW>
        cur = t;
    end
end
if strlength(cur) > 0, lines(end+1,1) = cur; end %#ok<AGROW>
s = strjoin(lines, newline);
end

function r = tilebox(ax, fig, pad)
if nargin<3, pad=[0 0 0 0]; end
px = getpixelposition(ax,true); pf = getpixelposition(fig,true);
r = [ (px(1)+pad(1))/pf(3), (px(2)+pad(2))/pf(4), ...
      max(0,(px(3)+pad(3))/pf(3)), max(0,(px(4)+pad(4))/pf(4)) ];
end

% ====== contrast helpers ======
function tz = zero_cross_after_peak(t, y)
% 主峰后首次与 0 的线性插值交点；找不到则 NaN
if all(~isfinite(y)), tz = NaN; return; end
[~, imx] = max(y);
tz = NaN;
k = imx + 1;
while k <= numel(t)
    if y(k-1) > 0 && y(k) <= 0
        % 线性插值
        x1 = t(k-1); x2 = t(k);
        y1 = y(k-1); y2 = y(k);
        if y2 == y1
            tz = x2;
        else
            tz = x1 + (0 - y1) * (x2 - x1) / (y2 - y1);
        end
        break
    end
    k = k + 1;
end
end

function tc = pos_centroid(t, y)
% 正瓣质心：∫ t·max(y,0) dt / ∫ max(y,0) dt
yp = max(y,0);
A  = trapz(t, yp);
if ~isfinite(A) || A <= eps
    tc = NaN;
else
    tc = trapz(t, yp .* t) / A;
end
end
