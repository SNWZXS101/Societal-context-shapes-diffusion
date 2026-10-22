function plot_panelB(out_root, features_mat, cluster_mat, varargin)
% PANEL B — 4x6 固定布局（4 原型 × 3 变异株 × {f,inc}）
% 每格独立轴：左列 f(t)，右列 inc(t)；空格自动做“同原型忽略变异株”的兜底代表。
%
% 示例：
% plot_panelB('results_cf','features.mat','smf_fpca_cluster_results.mat', ...
%   'Scenarios',{'NEG_STRONG','TGC_R2ZERO'}, ...
%   'VariantOrder',{'Alpha','Delta','Omicron'}, ...
%   'SavePNG',true,'OutName','panelB_4x6.png','FontSize',10);

% ---------- 参数 ----------
p = inputParser; p.FunctionName = mfilename;
p.addRequired('out_root',     @(x)ischar(x)||isstring(x));
p.addRequired('features_mat', @(x)ischar(x)||isstring(x));
p.addRequired('cluster_mat',  @(x)ischar(x)||isstring(x));
p.addParameter('Scenarios',   {'NEG_STRONG','TGC_R2ZERO'}, @(x)iscellstr(x)||all(string(x)==string(x)));
p.addParameter('KTarget',     3, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('VariantOrder',[], @(x)iscellstr(x)||isstring(x)||isempty(x));
p.addParameter('SavePNG',     false, @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('OutName',     'panelB_4x6.png', @(x)ischar(x)||isstring(x));
p.addParameter('FontSize',    10, @(x)isnumeric(x)&&isscalar(x)&&x>=6);
p.addParameter('LineWidth',   1.6, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('FigureSizePx',[1900 1200], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('FallbackIfEmpty', true, @(x)islogical(x)||ismember(x,[0 1]));
p.parse(out_root, features_mat, cluster_mat, varargin{:});
A = p.Results; out_root = char(out_root);

% ---------- 读数据 ----------
Sfeat = load(features_mat,'Features','f_cell');
Features = Sfeat.Features;  f_cell = Sfeat.f_cell;
N = height(Features); assert(numel(f_cell)==N, 'f_cell 与 Features 长度不一致');
Scf = load(fullfile(out_root,'cf_results.mat'),'Results'); Results = Scf.Results;

Sclu = load(cluster_mat);
[Rclu, labels, ~] = resolve_labels_and_proto_3x4(Sclu, A.KTarget, f_cell); % 带 u 的兜底
labels = labels(:); assert(numel(labels)==N, '聚类标签长度与 Features 不一致');

% 三列变异株
if isempty(A.VariantOrder)
    v_all = unique(string(Features.variant),'stable');
    variants = v_all(1:min(3,numel(v_all)));
else
    variants = string(A.VariantOrder(:));
end
if numel(variants)<3, variants(end+1:3) = "(empty)"; end
variants = variants(1:3);

% ---------- 代表波选择（4×3 + 兜底） ----------
Kshow = 4; Vshow = 3; scen_main = string(A.Scenarios(1));
rep_idx    = nan(Kshow,Vshow);
rep_idx_fb = nan(Kshow,1);   % 兜底：同原型不分变异株

have_proto = ~isempty(Rclu.proto) && size(Rclu.proto,1)==numel(Rclu.u);
% 先准备每个原型的“全体（不分变异株）”兜底代表
for k=1:Kshow
    idxs_k = find(labels==k);
    if isempty(idxs_k), continue; end
    target = get_target_curve(have_proto, Rclu, f_cell, idxs_k, k);
    [rep_idx_fb(k), ~] = pick_representative(out_root, Results, Features, f_cell, idxs_k, target, scen_main, Rclu.u);
end

% 再按 (k,v) 精确选择
for k=1:Kshow
  for v=1:Vshow
    if variants(v)=="(empty)"; continue; end
    idxs = find((labels==k) & strcmp(string(Features.variant), variants(v)));
    if isempty(idxs)
        if A.FallbackIfEmpty && isfinite(rep_idx_fb(k))
            rep_idx(k,v) = rep_idx_fb(k);   % 启用兜底
        end
        continue;
    end
    target = get_target_curve(have_proto, Rclu, f_cell, idxs, k);
    [rep_idx(k,v), ~] = pick_representative(out_root, Results, Features, f_cell, idxs, target, scen_main, Rclu.u);
  end
end

% ---------- 画布：4×6 ----------
f = figure('Color','w','Units','pixels','Position',[80 80 A.FigureSizePx]);
tlo = tiledlayout(4,6,'TileSpacing','compact','Padding','compact');
fs = A.FontSize;
% 颜色
Cbase  = [0 0 0];
Cproto = smf_colorlib('maina8');                    % 4 行颜色
if size(Cproto,1)<4, Cproto = repmat(Cproto,ceil(4/size(Cproto,1)),1); end
rowColor = @(k) Cproto(mod(k-1,size(Cproto,1))+1,:); % 第 k 行的颜色

% 场景线型（与之前一致：第1条实线，第2条虚线）
styles    = {'-','--',':','-.'};
pickStyle = @(s) styles{mod(s-1,numel(styles))+1};



% tile 索引：行 k，列 v：f 轴位于 col = 2*(v-1)+1；inc 轴位于 col+1
for k=1:Kshow
  for v=1:Vshow
    col_f = 2*(v-1)+1; col_i = col_f+1;
    % ---- f(t) 轴 ----
    axF = nexttile((k-1)*6 + col_f); hold(axF,'on'); box(axF,'on'); set(axF,'FontSize',fs);
    ii = rep_idx(k,v);

    title_suffix = '';
    if isnan(ii) || variants(v)=="(empty)"
        axis(axF,'off');
        title(axF, sprintf('C%d × %s (n=0)', k, char(variants(v))), ...
            'Interpreter','none','FontSize',fs,'Color',Cproto(k,:));
        % inc 轴也置空
        axI = nexttile((k-1)*6 + col_i); axis(axI,'off');
        continue;
    else
        % 如果是兜底来的，标注 fallback
        if A.FallbackIfEmpty
            mask_exact = (labels(ii)==k) & strcmp(string(Features.variant(ii)), variants(v));
            if ~mask_exact, title_suffix='  [fallback]'; end
        end
    end

    % 载入 TS
    wid = string(Features.wave_id(ii));
    TS  = load_TS_local(out_root, wid);

    % f(t)
    plot(axF, TS.t, TS.f, '-', 'Color', Cbase, 'LineWidth', 4);
    for s = 1:numel(A.Scenarios)
        tag = upper(string(A.Scenarios{s}));
        ff  = "f_"+tag;
        if isfield(TS, ff)
            plot(axF, TS.t, TS.(ff), pickStyle(s), ...
                'Color', rowColor(k), 'LineWidth', A.LineWidth+0.5);
        end
    end

    yl = ylim(axF); plot(axF,[TS.t(1) TS.t(end)],[0 0],'-','Color',[0 0 0 0.15],'LineWidth',1,'HandleVisibility','off'); ylim(axF,yl);
    xlabel(axF,'Day'); ylabel(axF,'f(t)');
    ttl = sprintf('%s-%s (%s | C%d)%s', string(Features.province(ii)), string(Features.wave_idx(ii)), ...
        char(variants(v)), k, title_suffix);
    title(axF, ttl, 'Interpreter','none','FontSize',fs,'Color',Cproto(k,:));

    % ---- inc(t) 轴 ----
    axI = nexttile((k-1)*6 + col_i); hold(axI,'on'); box(axI,'on'); set(axI,'FontSize',fs);
    plot(axI, TS.t, TS.inc, '-', 'Color', Cbase, 'LineWidth', 4);
    for s = 1:numel(A.Scenarios)
        tag = upper(string(A.Scenarios{s}));
        ff  = "inc_"+tag;
        if isfield(TS, ff)
            plot(axI, TS.t, TS.(ff), pickStyle(s), ...
                'Color', rowColor(k), 'LineWidth', A.LineWidth+0.5);
        end
    end

    xlabel(axI,'Day'); ylabel(axI,'inc');

    % %Avoided（主场景）
    pct_main = lookup_pct(Results, wid, scen_main);
    if isnan(pct_main), pct_main = compute_pct_from_TS(out_root, wid, char(scen_main)); end
    text(axI, 0.98, 0.95, sprintf('%.1f%%', pct_main), ...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top', ...
        'Color',[0.2 0.2 0.2],'FontSize',fs,'FontWeight','bold');
  end
end

% 顶部总图例
lgH = gobjects(0); lgT = strings(0,1);

% (a) 线型=场景（用灰色示例，避免跟行色冲突）
grayDemo = [0.25 0.25 0.25];
lgH(end+1) = plot(nan,nan,'-','Color',Cbase,'LineWidth',2.2);          lgT(end+1) = "Baseline";
for s = 1:min(2,numel(A.Scenarios))   % 只演示前两条
    lgH(end+1) = plot(nan,nan, pickStyle(s), 'Color', grayDemo, ...
        'LineWidth', A.LineWidth+0.5);
    lgT(end+1) = "Scenario: " + string(A.Scenarios{s});
end

% (b) 颜色=原型（用色块/方形标记示例）
for kk = 1:4
    lgH(end+1) = plot(nan,nan,'s','MarkerFaceColor',rowColor(kk), ...
        'MarkerEdgeColor',rowColor(kk), 'LineStyle','none','MarkerSize',8);
    lgT(end+1) = "Prototype C" + kk;
end

% 顶部放置
try
    lg = legend(tlo, lgH, lgT, 'Orientation','horizontal', 'NumColumns', 7, ...
        'Interpreter','none','Box','off');
    lg.Layout.Tile = 'north';
catch
    axAll = findobj(tlo,'Type','axes');
    if ~isempty(axAll)
        legend(axAll(1), lgH, lgT, 'Orientation','horizontal', 'NumColumns', 7, ...
            'Interpreter','none','Box','off','Location','northoutside');
    end
end


% 导出
if A.SavePNG
    out_png = fullfile(out_root,'figures',char(A.OutName));
    if ~isfolder(fileparts(out_png)), mkdir(fileparts(out_png)); end
    exportgraphics(gcf, out_png, 'Resolution', 300);
    fprintf('[plot_panelB] Saved → %s\n', out_png);
end
end

% ===== 代表选择/工具 =====
function target = get_target_curve(have_proto, Rclu, f_cell, idxs, k)
if have_proto && size(Rclu.proto,2)>=k
    target = Rclu.proto(:,k);
else
    Fc_tmp = resamp_many_to_u(f_cell(idxs), Rclu.u);
    target = median(Fc_tmp,2,'omitnan');
end
end

function [best_idx, scores] = pick_representative(out_root, Results, Features, f_cell, idxs, target, scen_main, u)
% 从 idxs 中选代表：形状距离(target)最小，且 %Avoided 接近组内中位
% 输入：
%   - out_root, Results, Features：读取 %Avoided 或从 TS 现算
%   - f_cell：原始 f(t) 列表
%   - idxs：候选样本的全局索引（对应 Features 的行）
%   - target：原型/中位曲线（在 u 上）
%   - scen_main：用于代表角标的主场景名（string）
%   - u：归一化时间网格（0–1）

% 1) 形状距离
Fc    = resamp_many_to_u(f_cell(idxs), u);                % M × m
dists = sum(abs(Fc - target), 1, 'omitnan')';             % m × 1  (L1)

% 2) %Avoided（主场景）
pcts = nan(numel(idxs),1);
for j = 1:numel(idxs)
    wid_j = string(Features.wave_id(idxs(j)));
    pct   = lookup_pct(Results, wid_j, scen_main);
    if isnan(pct)
        pct = compute_pct_from_TS(out_root, wid_j, char(scen_main));
    end
    pcts(j) = pct;
end
medp = median(pcts,'omitnan');

% 3) 排序：先形状最近，再 |pct-中位| 最小
[~, ord] = sortrows([dists, abs(pcts - medp)], [1 2]);
best_idx = idxs(ord(1));

if nargout > 1
    scores = struct('dist', dists(ord(1)), 'pct', pcts(ord(1)), 'med', medp);
end
end




function [R, labels, K] = resolve_labels_and_proto_3x4(S, Kreq, f_cell)
flds = fieldnames(S);
if numel(flds)==1 && isstruct(S.(flds{1})), S = S.(flds{1}); end
R = struct('u',[],'proto',[]);
labels=[]; K=[];

if isfield(S,'byK') && ~isempty(S.byK)
    KK = arrayfun(@(z) getfield_safe(z,'K',NaN), S.byK);
    pos = find(KK==Kreq,1);
    if ~isempty(pos)
        blk   = S.byK(pos);
        labels= pickf(blk, ["labels","assign","assignments","idx"]);
        R.u   = pickf(blk, ["u","grid","U"]);
        R.proto = pickf(blk, ["prototypes","proto","DBA","centers"]);
        if ~isempty(labels), K = max(labels); end
    end
end
if isempty(labels)
    clu     = pickf(S,'cluster');
    labels  = pickf(clu, ["labels","assign","assignments","idx"]);
    R.u     = pickf(S, ["u","grid","U"]); if isempty(R.u), R.u = pickf(clu, ["u","grid","U"]); end
    R.proto = pickf(clu, ["prototypes","proto","DBA","centers"]);
    if ~isempty(labels), K = max(labels); end
end
if ~isempty(R.proto)
    if size(R.proto,2)~=max(labels) && size(R.proto,1)==max(labels)
        R.proto = R.proto.';
    end
end
if isempty(R.u) || ~isvector(R.u)
    M = [];
    if ~isempty(R.proto), M = size(R.proto,1);
    elseif nargin>=3 && ~isempty(f_cell)
        L = cellfun(@numel, f_cell(:)); L = L(L>=2 & isfinite(L));
        if ~isempty(L), M = round(median(L)); end
    end
    if isempty(M) || M<2, M = 101; end
    R.u = linspace(0,1,M)'; 
else
    R.u = R.u(:);
end
end

function val = pickf(T, keys)
val = [];
if isempty(T), return; end
if ischar(keys) || isstring(keys), keys = string(keys); end
if istable(T)
    for i=1:numel(keys), k=char(keys(i));
        if any(strcmpi(T.Properties.VariableNames,k)), val=T.(k); return; end
    end
elseif isstruct(T)
    f = fieldnames(T);
    for i=1:numel(keys), k=char(keys(i));
        hit = f(strcmpi(f,k)); if ~isempty(hit), val=T.(hit{1}); return; end
    end
end
end

function v = getfield_safe(S, k, def), v = def; if isstruct(S)&&isfield(S,k)&&~isempty(S.(k)), v=S.(k); end, end

function Fc = resamp_many_to_u(f_list, u)
M=numel(u); m=numel(f_list); Fc=nan(M,m);
for j=1:m
    fj = double(f_list{j}(:)); Tj=numel(fj);
    if Tj<2, continue; end
    tj=(0:Tj-1)'; uj=(tj-tj(1))/(tj(end)-tj(1)+eps);
    Fc(:,j) = interp1(uj, fj, u(:), 'linear', 'extrap');
end
end

function TS = load_TS_local(out_root, wid)
p = fullfile(out_root,'timeseries', sprintf('wave_%s.mat', regexprep(char(wid),'[^a-zA-Z0-9-_]','_')));
if ~isfile(p), error('缺少 timeseries：%s', p); end
TSS = load(p,'TS'); TS = TSS.TS;
end

function pct = lookup_pct(Results, wid, tag)
pct = NaN;
mask = strcmpi(string(Results.wave_id), string(wid)) & strcmpi(string(Results.scenario), string(tag));
if any(mask)
    v = Results.avoided_pct(mask);
    if ~isempty(v), pct = double(v(1)); end
end
end

function pct = compute_pct_from_TS(out_root, wid, tag)
pct = NaN;
try
    TS = load_TS_local(out_root, wid);
    inc0 = TS.inc(:); incs = TS.("inc_"+upper(tag))(:);
    Cb = nansum(inc0); Cs = nansum(incs);
    pct = 100*(Cb - Cs)/max(Cb, eps);
catch
end
end
