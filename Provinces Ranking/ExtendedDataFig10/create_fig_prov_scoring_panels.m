function create_fig_prov_scoring_panels(scores_root, varargin)
% CREATE_FIG_PROV_SCORING_PANELS  (R2024b)
% Panels:
%   (A) RRI–MEI bubble matrix (area=SRS, color=Region; region from features_aug.mat)
%   (B) LPI waterfall (first→last gradient) with 95% CI from bootstrap
%   (C) Upward ridge plot of rank uncertainty (rows sorted by P(Top10))
%
% Upstream outputs (timestamp=false):
%   scores/map_ready_SRS.csv          % province, RRI_VA, MEI_VA, SRS, rank
%   scores/province_progress.csv      % province, first_SRS, last_SRS, slope, LPI, ...
%   scores/province_trajectories.csv  % province, order, SRS_wave, RRI_wave, MEI_wave
%
% Bootstrap (if enabled):
%   - uses province_trajectories.csv
%   - fixed scale mapping (min–max from baseline, consistent with compute_lpi)
%   - saves scores/lpi_boot.(mat|csv) with province names
%
% Example:
%   create_fig_prov_scoring_panels('scores', ...
%       'features_mat','features_aug.mat', ...
%       'do_bootstrap', true, 'B', 400, 'seed', 2025, ...
%       'slope_method','theilsen', ...   % 'theilsen' | 'robustfit'
%       'w_delta',0.5,'w_slope',0.5, ...
%       'export_pdf',false,'out_pdf','Fig_ScoringPanels.pdf');

% ---------- Args ----------
p = inputParser; p.FunctionName = mfilename;
addRequired(p,'scores_root',@(s)ischar(s)||isstring(s));
addParameter(p,'features_mat','features_aug.mat',@(s)ischar(s)||isstring(s));
addParameter(p,'do_bootstrap',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'B',400,@(x)isnumeric(x)&&isscalar(x)&&x>=50);
addParameter(p,'seed',2025,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'slope_method','theilsen',@(s)any(strcmpi(s,{'theilsen','robustfit'})));
addParameter(p,'w_delta',0.5,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'w_slope',0.5,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'export_pdf',false,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'out_pdf','Fig_ScoringPanels.pdf',@(s)ischar(s)||isstring(s));
parse(p,scores_root,varargin{:}); A = p.Results;

% ---------- Load core tables ----------
SRS = readtable(fullfile(A.scores_root,'map_ready_SRS.csv'),'VariableNamingRule','preserve');
LPI = readtable(fullfile(A.scores_root,'province_progress.csv'),'VariableNamingRule','preserve');
TRJ = readtable(fullfile(A.scores_root,'province_trajectories.csv'),'VariableNamingRule','preserve');
SRS.province = string(SRS.province);
LPI.province = string(LPI.province);
TRJ.province = string(TRJ.province);

% ---------- Province → Region from features_aug.mat ----------
Reg = build_region_from_features_local(unique(SRS.province,'stable'), A.features_mat);

% ---------- Bootstrap (fixed scale; aligned with compute_lpi) ----------
Boot = struct();
if A.do_bootstrap
    Bargs = struct('outdir',A.scores_root,'B',A.B,'seed',A.seed, ...
                   'slope_method',A.slope_method, 'w_delta',A.w_delta,'w_slope',A.w_slope);
    Boot = bootstrap_lpi_and_ranks(Bargs);            % fixed scale + time order
    save_lpi_boot(Boot, A.scores_root);               % includes province names
end

% ---------- Figure canvas ----------
fig = figure('Color','w','Position',[80 80 1400 900],'WindowStyle','normal');
movegui(fig,'center'); plotedit(fig,'on');  % allow drag/resize
axA = axes('Parent',fig,'Units','normalized','Position',[0.08 0.55 0.58 0.40]); % Panel A
axB = axes('Parent',fig,'Units','normalized','Position',[0.70 0.55 0.26 0.40]); % Panel B
axC = axes('Parent',fig,'Units','normalized','Position',[0.08 0.08 0.88 0.38]); % Panel C

% ---------- Draw panels ----------
draw_panelA_bubble_matrix(axA, innerjoin(SRS,Reg,'Keys','province'));
draw_panelB_lpi_waterfall(axB, LPI, Boot);
draw_panelC_rank_ridge(axC, Boot);

% ---------- Export ----------
if A.export_pdf
    set(fig,'Renderer','painters');
    print(fig, A.out_pdf, '-dpdf','-bestfit');
    fprintf('[export] saved: %s\n', A.out_pdf);
end
end

% ===================== Panel A =====================
function draw_panelA_bubble_matrix(ax, T)
axes(ax); cla(ax); hold(ax,'on');

% Region colors (7 high-contrast "sunrise → royal" palette)
regList = categories(categorical(T.region));
K = numel(regList);
C7 = [0.98 0.78 0.15; 0.96 0.62 0.00; 0.89 0.43 0.32; ...
      0.78 0.33 0.54; 0.42 0.36 0.90; 0.25 0.47 0.71; 0.15 0.42 0.55];
if K > 7
    % extend if >7
    Ck = extend_palette_rows(C7, K);
else
    Ck = C7(1:K,:);
end
cmap = containers.Map(cellstr(regList), num2cell(Ck,2));

% Bubble size from SRS (area in pt^2)
srs = double(T.SRS); srs = max(srs,0);
Amin=50; Amax=700;
Aarea = Amin + (srs - min(srs))/max(1,(max(srs)-min(srs))) * (Amax - Amin);

% Scatter
n = height(T); col = zeros(n,3);
for i=1:n, col(i,:)=cmap(char(T.region(i))); end
scatter(ax, T.RRI_VA, T.MEI_VA, Aarea, col, 'filled', ...
        'MarkerEdgeColor',[0.18 0.18 0.18], 'LineWidth',0.6, 'MarkerFaceAlpha',0.9);

% Quadrant reference
xline(ax,50,':','Linewidth',2.5,'Color',[0.55 0.55 0.55]); yline(ax,50,':','Linewidth',2.5,'Color',[0.55 0.55 0.55]);

% Labels & axes
xlabel(ax,'RRI_{VA}'); ylabel(ax,'MEI_{VA}');
title(ax,'(A) Provincial RRI–MEI bubble matrix (area=SRS; color=Region)','FontSize',13);
xlim(ax,[0 100]); ylim(ax,[0 100]); grid(ax,'on'); set(ax,'Layer','top'); box(ax,'off'); set(ax,'TickLength',[0 0]);

% Abbreviations under bubbles
abbr = abbrev_province(T.province);
for i=1:n
    text(ax, T.RRI_VA(i), T.MEI_VA(i)-2.3, abbr(i), 'HorizontalAlignment','center', ...
        'VerticalAlignment','top','FontSize',10,'Color',[0.20 0.20 0.20]);
end

% Legend
uRegs = unique(T.region,'stable'); ph = gobjects(numel(uRegs),1);
for k=1:numel(uRegs)
    ph(k) = plot(ax, nan,nan,'o','MarkerFaceColor',cmap(char(uRegs(k))), ...
        'MarkerEdgeColor',[0.18 0.18 0.18],'LineWidth',0.6,'MarkerSize',10);
end
legend(ax, ph, cellstr(uRegs), 'Location','northwest','Box','off');

% Right-side annotation: abbreviation ↔ full name
fig = ancestor(ax,'figure');
lst = strings(n,1);
for i=1:n, lst(i) = sprintf('%-8s — %s', abbr(i), string(T.province(i))); end
lst = join(lst, newline);
annotation(fig,'textbox','Units','normalized', ...
    'Position',[0.86 0.12 0.12 0.78], 'String',lst, 'FontSize',10, ...
    'EdgeColor',[0.85 0.85 0.85], 'BackgroundColor',[1 1 1], ...
    'Interpreter','none', 'Margin',6);
end

% ===================== Panel B =====================
function draw_panelB_lpi_waterfall(ax, LPI, Boot)
% Horizontal waterfall: first→last gradient, black dot at LPI, 95% CI from bootstrap.
axes(ax); cla(ax); hold(ax,'on');

LPI.province = string(LPI.province);
[~,ord] = sort(LPI.LPI, 'descend', 'MissingPlacement','last');
LPI = LPI(ord,:); 
prov = string(LPI.province);
N = height(LPI);

% Align CI from Boot if present
ciL = nan(N,1); ciU = ciL;
if isstruct(Boot) && isfield(Boot,'lpi_ci') && isfield(Boot,'province_order')
    provB = string(Boot.province_order);
    [tf,loc] = ismember(prov, provB);
    ciL(tf) = Boot.lpi_ci(loc(tf),1);
    ciU(tf) = Boot.lpi_ci(loc(tf),2);
end

% Draw horizontal gradient strips (first→last) + black dot at LPI + CI whiskers
y = 1:N;
for i=1:N
    y0 = y(i);

    % gradient strip: first SRS (light) -> last SRS (dark)
    x1 = LPI.first_SRS(i); 
    x2 = LPI.last_SRS(i);
    X = [x1 x2 x2 x1]; 
    Y = [y0-0.36 y0-0.36 y0+0.36 y0+0.36];
    c_first = [0.88 0.90 0.98];   % light
    c_last  = [0.20 0.40 0.90];   % darker
    fvc = [c_first; c_last; c_last; c_first];
    patch('XData',X,'YData',Y,'FaceVertexCData',fvc, 'FaceColor','interp', ...
          'EdgeColor',[0.25 0.25 0.25],'LineWidth',0.7);

    % 95% CI whiskers
    if isfinite(ciL(i)) && isfinite(ciU(i))
        line([ciL(i) ciU(i)], [y0 y0], 'Color',[0.15 0.15 0.15], 'LineWidth',1.2);
        plot([ciL(i) ciL(i)], [y0-0.18 y0+0.18], '-', 'Color',[0.15 0.15 0.15], 'LineWidth',1.2);
        plot([ciU(i) ciU(i)], [y0-0.18 y0+0.18], '-', 'Color',[0.15 0.15 0.15], 'LineWidth',1.2);
    end

    % LPI point estimate: uniform black dot
    plot(LPI.LPI(i), y0, 'ko', 'MarkerFaceColor','k', 'MarkerSize',5.5, 'LineWidth',0.5);

    % numeric label
    text(LPI.LPI(i)+1.2, y0+0.5, sprintf('%.1f', LPI.LPI(i)), ...
         'HorizontalAlignment','left','VerticalAlignment','middle', ...
         'Color',[0.25 0.25 0.25],'FontSize',10.5);
end

% Axes styling
set(ax,'YTick',y,'YTickLabel',prov,'YDir','reverse','TickLength',[0 0], 'Layer','top','Box','off');
xlabel(ax,'LPI'); 
ylim(ax,[0.5 N+0.5]);

xmin = max(0, min([0;ciL;LPI.LPI])-3);
xmax = min(100, max([100;ciU;LPI.LPI])+6);
xlim(ax,[xmin xmax]);

% Title & grid
title(ax, {'LPI waterfall with 95% CI', ...
           ['Black dot = point estimate; ' ...
           'strip gradient = First SRS \rightarrow Last SRS']});
grid(ax,'on');

% Mini legend for the strip gradient (NOT a colorbar)
% --- Mini legend for the strip gradient (NOT a colorbar) — place at top-left ---
xL = ax.XLim(1); xR = ax.XLim(2);
yB = ax.YLim(2); yT = ax.YLim(1);

% 边距（相对轴域尺寸）
dx = 0.05*(xR - xL);
dy = 0.1*(yT - yB);

% 左上角起点
x0 = xL + dx;
y0 = yT - dy;

% 色块尺寸（数据坐标单位，可按需要微调）
w = 4; 
h = 1;

c_first = [0.88 0.90 0.98];   % 与条带左端一致（First SRS）
c_last  = [0.20 0.40 0.90];   % 与条带右端一致（Last SRS）

% 第一行：First SRS（靠上）
patch([x0 x0+w x0+w x0], [y0-h y0-h y0 y0], c_first, ...
      'EdgeColor',[.25 .25 .25], 'LineWidth',0.6);
text(x0+w+0.6, y0-0.5*h, 'First SRS', ...
     'HorizontalAlignment','left','VerticalAlignment','middle', ...
     'FontSize',11, 'Color',[.25 .25 .25]);

% 第二行：Last SRS（在其下方）
y1 = y0 - (h + 0.35*h);  % 行距可调（这里取 0.35*h）
patch([x0 x0+w x0+w x0], [y1-h y1-h y1 y1], c_last, ...
      'EdgeColor',[.25 .25 .25], 'LineWidth',0.6);
text(x0+w+0.6, y1-0.5*h, 'Last SRS', ...
     'HorizontalAlignment','left','VerticalAlignment','middle', ...
     'FontSize',11, 'Color',[.25 .25 .25]);


end



% ===================== Panel C =====================
function draw_panelC_rank_ridge(ax, Boot)
% Upward ridges; sorted by P(Top10); row tint = sunrise→slate (no colorbar).
axes(ax); cla(ax); hold(ax,'on');

if ~(isstruct(Boot) && isfield(Boot,'rank_samples') && isfield(Boot,'province_order'))
    text(0.5,0.5,'Bootstrap disabled (no rank_samples)','Units','normalized','HorizontalAlignment','center');
    title(ax,'(C) Upward ridge plot of rank uncertainty (sorted by P(Top10))');
    axis(ax,'off'); return;
end

R = Boot.rank_samples;              % P×B
prov = string(Boot.province_order); % P×1
P = numel(prov);

% 排序：按 P(Top10) 从大到小
pt10 = mean(R<=10, 2, 'omitnan');
[~, ord] = sort(pt10, 'descend', 'MissingPlacement','last');
prov = prov(ord);
R    = R(ord,:);
pt10 = pt10(ord);

% 行配色（仅编码行序）
Crows = ridge_row_cmap_local(P);

% —— 关键：坐标轴向上增长，YTick 递增 —— 
set(ax,'YDir','normal');   
ystep   = 1.08;
yposAsc = (1:P) * ystep;            % 递增 tick 位置（底→顶）

% 画 ridges：把第 i 条画在 y = (P-i+1)*ystep（顶端最大）
xgrid = linspace(1, P, 600);
for i = 1:P
    ri = R(i,:); ri = ri(isfinite(ri));
    if numel(ri) < 3, continue; end
    [f, xi] = ksdensity(ri, xgrid, 'Bandwidth', 0.6);
    f = f / max(f);
    y0 = (P - i + 1) * ystep;       % 越前的省越靠上
    patch([xi, fliplr(xi)], [y0 + f, fliplr(y0 + zeros(size(f)))], ...
          Crows(i,:), 'EdgeColor', darken_local(Crows(i,:),0.25), ...
          'FaceAlpha', 0.82, 'LineWidth', 0.8);
end

% 轴与标注：YTick 递增，但标签反转以使最上面显示第1名
set(ax, 'YTick', yposAsc, 'YTickLabel', flipud(prov), ...
        'TickLength',[0 0], 'Box','off', 'Layer','top');

% 在右侧加 P(Top10) 数值列
pt10_flip = flipud(pt10);             % 与 YTickLabel 对齐
xRightSep = P + 0.3;                  % 分隔线位置
xRightTxt = P + 0.6;                  % 数值文本位置
xlim(ax,[1 P+1.8]);                   % 预留右侧空间

xline(ax, xRightSep, '-', 'Color', [0.8 0.8 0.85], 'LineWidth', 1.0);
for j = 1:P
    text(ax, xRightTxt, yposAsc(j), sprintf('%.0f%%', 100*pt10_flip(j)), ...
        'HorizontalAlignment','left','VerticalAlignment','middle', ...
        'Color',[0.25 0.25 0.25], 'FontSize', 11);
end
% 右上角列名
text(ax, xRightTxt, yposAsc(end)+0.6*ystep, 'P(Top10)', ...
    'HorizontalAlignment','left','VerticalAlignment','bottom', ...
    'FontWeight','bold','Color',[0.20 0.20 0.20]);

xlabel(ax,'Rank (1 = best)');
title(ax,'(C) Upward ridge plot of rank uncertainty (sorted by P(Top10))');
grid(ax,'on'); set(ax,'GridAlpha',0.08);
xline(ax, 10.5, ':', 'Color', [0.45 0.45 0.45]); % Top-10 阈值
ylim(ax,[ystep*0.3, ystep*(P+0.9)]);
end




% ===================== Bootstrap (fixed scale) =====================
function OUT = bootstrap_lpi_and_ranks(A)
% Uses scores/province_trajectories.csv; fixed scale = baseline min–max across provinces.
if ~isfield(A,'outdir'),  A.outdir='scores'; end
if ~isfield(A,'B'),       A.B=400;           end
if ~isfield(A,'seed'),    A.seed=2025;       end
if ~isfield(A,'w_delta'), A.w_delta=0.5;     end
if ~isfield(A,'w_slope'), A.w_slope=0.5;     end
if ~isfield(A,'slope_method'), A.slope_method='theilsen'; end
rng(A.seed,'twister');

T = readtable(fullfile(A.outdir,'province_trajectories.csv'),'VariableNamingRule','preserve');
T.province = string(T.province);
need = {'province','order','SRS_wave'};
miss = setdiff(need, T.Properties.VariableNames);
if ~isempty(miss), error('[bootstrap] province_trajectories missing: %s', strjoin(miss,', ')); end

% numeric order
if ~isnumeric(T.order)
    ordNum = str2double(string(T.order));
    if any(isnan(ordNum))
        T = sortrows(T,'province');
        T.order = within_group_seq(T.province);
    else
        T.order = ordNum;
    end
end
T = sortrows(T,{'province','order'});

prov = unique(T.province,'stable'); P = numel(prov); B = A.B;

% ---- Baseline delta/slope (fixed bounds) ----
delta0 = nan(P,1); slope0 = nan(P,1);
for i=1:P
    Ti = T(T.province==prov(i),:); [~,ix]=sort(Ti.order); y = Ti.SRS_wave(ix);
    if numel(y)<2, continue; end
    delta0(i) = delta_ceiling_norm_local(y);
    if strcmpi(A.slope_method,'robustfit')
        b = robustfit((1:numel(y))', y); slope0(i) = b(2);
    else
        slope0(i) = theilsen_local((1:numel(y))', y);
    end
end
dmin = min(delta0,[],'omitnan'); dmax = max(delta0,[],'omitnan');
smin = min(slope0,[],'omitnan'); smax = max(slope0,[],'omitnan');

% ---- Bootstrap with time-order resampling ----
lpi_boot = nan(P,B); rank_samples = nan(P,B);
for b=1:B
    d_all = nan(P,1); s_all = nan(P,1);
    for i=1:P
        Ti = T(T.province==prov(i),:); [~,ix]=sort(Ti.order); y = Ti.SRS_wave(ix);
        n = numel(y); if n<2, continue; end
        idx = sort(randi(n,n,1)); yb = y(idx);  % preserve time order
        d_all(i) = delta_ceiling_norm_local(yb);
        if strcmpi(A.slope_method,'robustfit')
            bb = robustfit((1:numel(yb))', yb); s_all(i) = bb(2);
        else
            s_all(i) = theilsen_local((1:numel(yb))', yb);
        end
    end
    comp_delta = map_to_0_100_fixbounds(d_all, dmin, dmax);
    comp_slope = map_to_0_100_fixbounds(s_all, smin, smax);
    lpi_vec = A.w_delta*comp_delta + A.w_slope*comp_slope;
    lpi_boot(:,b) = lpi_vec;
    rank_samples(:,b) = tiedrank(-lpi_vec);
end

ciL = prctile(lpi_boot, 2.5, 2);
ciU = prctile(lpi_boot,97.5, 2);
OUT = struct('lpi_boot',lpi_boot,'province_order',prov(:), ...
             'lpi_ci',[ciL ciU], 'rank_samples',rank_samples);
end

function save_lpi_boot(Boot, outdir)
if ~isfolder(outdir), mkdir(outdir); end
X = Boot.lpi_boot; prov = string(Boot.province_order); [P,B] = size(X);
% .mat
save(fullfile(outdir,'lpi_boot.mat'), 'X','prov','B');
% .csv
lab = compose('sample_%03d',1:B);
T = [table(prov,'VariableNames',{'province'}), array2table(X,'VariableNames',lab)];
writetable(T, fullfile(outdir,'lpi_boot.csv'));
fprintf('[save_lpi_boot] Saved %s | %s  (P=%d, B=%d)\n', ...
    fullfile(outdir,'lpi_boot.mat'), fullfile(outdir,'lpi_boot.csv'), P, B);
end

% ===================== Region / Colors / Text helpers =====================
function Reg = build_region_from_features_local(provNeeded, features_mat)
S = load(features_mat); cand={'Features','Tfeat','features','T','Features_aug'}; hit='';
for nm=cand, if isfield(S,nm{1}), hit=nm{1}; break; end, end
assert(~isempty(hit),'features_mat: main table not found.');
T = S.(hit); assert(istable(T),'features_mat main object is not a table.');
V = string(T.Properties.VariableNames); L = lower(V);
col_p = find_col(V,L,{'^province$','^prov(ince)?$'},'province');
col_r = find_col(V,L,{'^region$','^geo[_]?region$','^macro[_]?region$','^area$'},'region');
P = string(T.(col_p)); R = string(T.(col_r));
[~,ia]=unique([P R],'rows','stable');
Map = table(P(ia), R(ia), 'VariableNames',{'province','region'});
need = string(unique(provNeeded)); [tf,loc]=ismember(need, Map.province);
assert(all(tf), 'features_aug region mapping missing provinces: %s', strjoin(cellstr(need(~tf)),', '));
Reg = Map(loc,:); ureg = unique(Reg.region); fprintf('[region] from features_mat: %d regions.\n', numel(ureg));
end
function real = find_col(V,L,pat,fallback)
real=""; for ii=1:numel(pat)
    m = ~cellfun('isempty', regexp(L, pat{ii}, 'once'));
    if any(m), real = V(find(m,1,'first')); return; end
end
if ~isempty(fallback), m = contains(L,fallback); if any(m), real=V(find(m,1,'first')); end, end
end

function ab = abbrev_province(prov)
% Pinyin-initial-like dotted abbreviations (adjust/extend map as needed)
M = containers.Map;
M('Beijing')='B.J.'; M('Tianjin')='T.J.'; M('Hebei')='H.B.'; M('Shanxi')='S.X.'; M('Neimenggu')='N.M.G.';
M('Liaoning')='L.N.'; M('Jilin')='J.L.'; M('Heilongjiang')='H.L.J.';
M('Shanghai')='S.H.'; M('Jiangsu')='J.S.'; M('Zhejiang')='Z.J.'; M('Anhui')='A.H.';
M('Fujian')='F.J.'; M('Jiangxi')='J.X.'; M('Shandong')='S.D.';
M('Henan')='H.N.'; M('Hubei')='H.Bei.'; M('Hunan')='H.Nan.';
M('Guangdong')='G.D.'; M('Guangxi')='G.X.'; M('Hainan')='H.N.';
M('Chongqing')='C.Q.'; M('Sichuan')='S.C.'; M('Guizhou')='G.Z.'; M('Yunnan')='Y.N.';
M('Shaanxi')='S.Xi.'; M('Gansu')='G.S.'; M('Qinghai')='Q.H.'; M('Ningxia')='N.X.'; 
M('Xinjiang')='X.J.'; M('Xizang')='X.Z.';
ab = strings(numel(prov),1);
for i=1:numel(prov)
    k = char(prov(i));
    if isKey(M,k), ab(i)=M(k); else, ab(i)=fallback_abbrev(k); end
end
end
function s = fallback_abbrev(name)
tok = regexp(name,'[A-Z][a-z]*','match'); if isempty(tok), tok={name}; end
ini = cellfun(@(x)x(1), tok); s = join(string(upper(ini(1:min(3,numel(ini))))),'.') + ".";
end

function C = ridge_row_cmap(N)
% "sunrise → slate" gradient for rows
a=[0.98 0.85 0.30]; b=[0.30 0.45 0.70]; t=linspace(0,1,N)'; C=a+(b-a).*t;
end
function y = darken(c,d), y=(1-d).*c; end

function seq = within_group_seq(g)
g = string(g); seq = zeros(size(g)); u=unique(g,'stable');
for i=1:numel(u), idx=find(g==u(i)); seq(idx)=1:numel(idx); end
end

% ===== compute_lpi-aligned helpers =====
function d = delta_ceiling_norm_local(y)
y = y(:); y = y(~isnan(y));
if numel(y)<2, d=NaN; return; end
first=y(1); last=y(end); delta=last-first;
if delta>=0, cap=max(1,100-first); d=100*delta/cap;
else,        cap=max(1,first);     d=100*delta/cap; end
end
function s = theilsen_local(x,y)
x=x(:); y=y(:); m=isfinite(x)&isfinite(y); x=x(m); y=y(m);
n=numel(x); if n<2, s=NaN; return; end
S=[]; 
for i=1:n-1
    dx=x(i+1:end)-x(i); dy=y(i+1:end)-y(i); m2=dx~=0;
    if any(m2), S=[S; dy(m2)./dx(m2)]; end %#ok<AGROW>
end
if isempty(S), s=NaN; else, s=median(S,'omitnan'); end
end
function z = map_to_0_100_fixbounds(x, xmin, xmax)
z = 50 + 0.*x;
if ~isfinite(xmin) || ~isfinite(xmax) || abs(xmax-xmin)<1e-12, return; end
z = 100 * (x - xmin) / (xmax - xmin);
z = min(max(z,0),100);
end

function Ck = extend_palette_rows(C, K)
% repeat base rows with slight lightness/saturation alternation
n=size(C,1); Ck=zeros(K,3);
for i=1:K
    base=C(mod(i-1,n)+1,:); cyc=floor((i-1)/n);
    if mod(cyc,2)==0, Ck(i,:)=min(1, base+0.08); else, Ck(i,:)=max(0, base-0.08); end
end
end


function C = ridge_row_cmap_local(N)
% Row tint for ridge plot: sunrise (gold) -> slate blue
a = [0.98 0.85 0.30];
b = [0.30 0.45 0.70];
t = linspace(0,1,N)';
C = a + (b - a) .* t;
end

function y = darken_local(c, d)
% Slightly darken an RGB color (0<=d<1)
y = (1 - d) .* c;
end



function C = colormap_factory_local(style, N)
% Discrete single-hue colormap, light->dark (top->bottom aesthetics).
if nargin<2, N=30; end
style = lower(string(style));
switch style
    case "orange", baseH = 28/360;  s0=0.25; s1=0.95; v0=1.00; v1=0.55;
    case "purple", baseH = 270/360; s0=0.20; s1=0.90; v0=0.98; v1=0.55;
    case "blue",   baseH = 210/360; s0=0.05; s1=0.85; v0=0.98; v1=0.55;
    case "green",  baseH = 130/360; s0=0.20; s1=0.90; v0=0.98; v1=0.55;
    case "red",    baseH = 5/360;   s0=0.25; s1=0.95; v0=0.98; v1=0.55;
    case "yellow", baseH = 48/360;  s0=0.05; s1=0.80; v0=1.00; v1=0.85;
    case "teal",   baseH = 170/360; s0=0.05; s1=0.85; v0=0.98; v1=0.60;
    otherwise, error('Unknown style: %s', style);
end
s = linspace(s0, s1, N)'; 
v = linspace(v0, v1, N)'; 
h = baseH * ones(N,1);
C = hsv2rgb([h, s, v]);
end
