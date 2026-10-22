function plot_panelD_heat(out_root, varargin)
% PLOT_PANELD_HEAT — Panel D 热图（NEG 与 TGC(t0 shift)），上下两块
% 数据口径与 counterfactual_main_heat 配套：
%   - Heat.NEG / Heat.TGC_T0 为长表，含列：wave_id, (eta|dt0), avoided_pct, cluster, ...
%   - Heat.Meta 仅含 variant/cluster 向量；此处不再从 Meta 取 cluster。
%
% 用法示例（推荐）：
%   plot_panelD_heat('results_cf_heat', ...
%       'GroupBy','variant_cluster', 'FeaturesMat','features.mat', ...
%       'Colormap','RdBu', 'SharedCLim',true, ...
%       'Annotate',true, 'AnnotateQuantile',0.90, 'AnnotateTopK',60, ...
%       'SavePNG',true,'OutName','panelD_heat');

% ---------- 参数 ----------
p = inputParser; p.FunctionName = mfilename;
addRequired(p,'out_root',@(x)ischar(x)||isstring(x));
addParameter(p,'GroupBy','variant',@(s) any(strcmpi(s,{'variant','cluster','variant_cluster'})));
addParameter(p,'FeaturesMat','',@(x)ischar(x)||isstring(x));  % 用于补齐 variant 映射（必须含 wave_id, variant）
addParameter(p,'Annotate',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'AnnotateFmt','%.1f',@(x)ischar(x)||isstring(x));
addParameter(p,'AnnotateTopK',0,@(x)isnumeric(x)&&isscalar(x)&&x>=0);           % 全局最多标注数
addParameter(p,'AnnotateQuantile',1.0,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<=1);% 仅标注 |v|≥分位数
addParameter(p,'FontSize',10,@(x)isnumeric(x)&&isscalar(x)&&x>=6);
addParameter(p,'SavePNG',false,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'SavePDF',false,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'OutName','panelD_heat',@(x)ischar(x)||isstring(x)); % 不带扩展名
addParameter(p,'Colormap','RdBu',@(x)ischar(x)||isstring(x));       % 默认用更“论文范”的 RdBu
addParameter(p,'SharedCLim',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'ClusterOrder','ascend', ...   % 'ascend'|'descend'|numeric vector
    @(x)ischar(x)||isstring(x)||isnumeric(x));


parse(p,out_root,varargin{:});
A = p.Results;

out_root = char(out_root);
heat_path = fullfile(out_root,'heat_results.mat');
if ~isfile(heat_path)
    error('找不到结果文件：%s（请先运行 counterfactual_main_heat）', heat_path);
end

S = load(heat_path,'Heat');
if ~isfield(S,'Heat'), error('heat_results.mat 中缺少 Heat 结构。'); end
Heat = S.Heat;

% ---------- 参数向量 ----------
if ~isfield(Heat,'Param') || ~all(isfield(Heat.Param, {'eta_list','dt0_list'}))
    error('Heat.Param 中缺少 eta_list / dt0_list。');
end
eta_list = unique(double(Heat.Param.eta_list(:)),'stable');     % NEG 横轴
dt0_list = unique(double(Heat.Param.dt0_list(:)),'stable');     % TGC 横轴

% ---------- 表 ----------
assert(isfield(Heat,'NEG') && istable(Heat.NEG), 'Heat.NEG 缺失或不是表');
assert(isfield(Heat,'TGC_T0') && istable(Heat.TGC_T0), 'Heat.TGC_T0 缺失或不是表');
Tneg = Heat.NEG;    % 期望列：wave_id, eta, avoided_pct, cluster
Ttgc = Heat.TGC_T0; % 期望列：wave_id, dt0, avoided_pct, cluster
need_neg = {'wave_id','eta','avoided_pct','cluster'};
need_tgc = {'wave_id','dt0','avoided_pct','cluster'};
assert(all(ismember(need_neg, Tneg.Properties.VariableNames)), 'Heat.NEG 缺列（需要 wave_id, eta, avoided_pct, cluster）');
assert(all(ismember(need_tgc, Ttgc.Properties.VariableNames)), 'Heat.TGC_T0 缺列（需要 wave_id, dt0, avoided_pct, cluster）');



% ---------- 变体映射（variant）来自 FeaturesMat ----------
variant_map = struct('ok',false);
if ~isempty(A.FeaturesMat)
    fpath = char(A.FeaturesMat);
    if isfile(fpath)
        Sf = load(fpath,'Features');
        if isfield(Sf,'Features') && istable(Sf.Features) ...
                && ismember('wave_id', Sf.Features.Properties.VariableNames)
            vcol = detect_alt_col(Sf.Features, "variant");
            if ~isempty(vcol)
                VarMap = Sf.Features(:, {'wave_id', vcol});
                VarMap.wave_id = string(VarMap.wave_id);
                VarMap.Properties.VariableNames{2} = 'variant';
                variant_map.ok = true;
                variant_map.wave = VarMap.wave_id;
                variant_map.var  = string(VarMap.variant);
            end
        end
    end
end
if ~variant_map.ok && any(strcmpi(A.GroupBy, {'variant','variant_cluster'}))
    warning('未提供有效的 FeaturesMat（含 wave_id 与 variant）。GroupBy=%s 将降级处理。', char(A.GroupBy));
end


% ---------- 为两个表各自构造分组键 ----------
[g_neg, g_tgc, grp_labels, ytitle] = build_group_keys(Tneg, Ttgc, A.GroupBy, variant_map);

% >>> 在这里按聚类重排行标签（若 GroupBy 涉及 cluster）
grp_labels = apply_cluster_order(grp_labels, A.ClusterOrder, A.GroupBy);




% ---------- 组×参数 的中位矩阵 ----------
Mneg = build_median_matrix_core(Tneg, 'eta',  eta_list, g_neg, grp_labels);
Mtgc = build_median_matrix_core(Ttgc, 'dt0', dt0_list, g_tgc, grp_labels);

% ---------- 画图（上下布局） ----------
f = figure('Color','w','Position',[120 80 1100 820], 'Renderer','opengl');
tlo = tiledlayout(f,2,1,'TileSpacing','compact','Padding','compact');

% 颜色 & caxis
cmap = get_nature_cmap(A.Colormap, 256);
if A.SharedCLim
    vmin = min([Mneg(:); Mtgc(:)], [], 'omitnan');
    vmax = max([Mneg(:); Mtgc(:)], [], 'omitnan');
    % 若含正负，零居中对称
    if isfinite(vmin) && isfinite(vmax) && vmin < 0 && vmax > 0
        m = max(abs([vmin, vmax])); vmin = -m; vmax = m;
    end
    if ~isfinite(vmin)||~isfinite(vmax) || vmin==vmax, vmin=0; vmax=100; end
else
    vmin=[]; vmax=[];
end

% ---- 上：NEG ----
ax1 = nexttile(tlo,1); hold(ax1,'on');
imagesc(ax1, eta_list, 1:numel(grp_labels), Mneg);
set(ax1,'YDir','normal','FontSize',A.FontSize);
set(ax1,'YTick',1:numel(grp_labels),'YTickLabel',cellstr(grp_labels));
xlabel(ax1,'\eta (negative-phase amplification)');
ylabel(ax1, ytitle);
title(ax1,'%Avoided — NEG');
colormap(ax1,cmap);
if ~isempty(vmin), caxis(ax1,[vmin vmax]); end
himg1 = findobj(ax1,'Type','Image'); set(himg1,'AlphaData',~isnan(Mneg));
grid(ax1,'on'); ax1.Layer='top'; ax1.LineWidth=0.8;
cb1 = colorbar(ax1,'Location','eastoutside'); cb1.Label.String = '% avoided'; cb1.Label.FontSize = A.FontSize;

if A.Annotate
    local_annotate_heat(ax1, Mneg, eta_list, 1:numel(grp_labels), A.AnnotateFmt, max(8,A.FontSize-1), A.AnnotateTopK, A.AnnotateQuantile);
end

% ---- 下：TGC (t0 shift) ----
ax2 = nexttile(tlo,2); hold(ax2,'on');
imagesc(ax2, dt0_list, 1:numel(grp_labels), Mtgc);
set(ax2,'YDir','normal','FontSize',A.FontSize);
set(ax2,'YTick',1:numel(grp_labels),'YTickLabel',cellstr(grp_labels));
xlabel(ax2,'\Delta t_0 (days)');
ylabel(ax2, ytitle);
title(ax2,'%Avoided — TGC (t_0 shift)');
colormap(ax2,cmap);
if ~isempty(vmin), caxis(ax2,[vmin vmax]); end
himg2 = findobj(ax2,'Type','Image'); set(himg2,'AlphaData',~isnan(Mtgc));
grid(ax2,'on'); ax2.Layer='top'; ax2.LineWidth=0.8;
cb2 = colorbar(ax2,'Location','eastoutside'); cb2.Label.String = '% avoided'; cb2.Label.FontSize = A.FontSize;

if A.Annotate
    local_annotate_heat(ax2, Mtgc, dt0_list, 1:numel(grp_labels), A.AnnotateFmt, max(8,A.FontSize-1), A.AnnotateTopK, A.AnnotateQuantile);
end

% ---- 导出 ----
base = char(A.OutName);
if A.SavePNG
    out_png = fullfile(out_root, base + ".png");
    try
        exportgraphics(f, out_png, 'Resolution', 300, 'BackgroundColor','white');
    catch
        print(f, out_png, '-dpng', '-r300');
    end
    fprintf('[plot_panelD_heat] Saved PNG → %s\n', out_png);
end
if A.SavePDF
    out_pdf = fullfile(out_root, base + ".pdf");
    try
        exportgraphics(f, out_pdf, 'ContentType','image', 'BackgroundColor','white', 'Resolution', 300);
    catch
        print(f, out_pdf, '-dpdf', '-r300');
    end
    fprintf('[plot_panelD_heat] Saved PDF → %s\n', out_pdf);
end
end

% ========== 从两个表构造分组键（关键修复：cluster 取自表；variant 由 FeaturesMat 映射） ==========
function [g_neg, g_tgc, grp_labels, ytitle] = build_group_keys(Tneg, Ttgc, groupBy, vmap)
gb = lower(string(groupBy));
switch gb
    case "cluster"
        g_neg = label_cluster(Tneg.cluster);
        g_tgc = label_cluster(Ttgc.cluster);
        grp_labels = union_unique(g_neg, g_tgc);
        ytitle = "CLUSTER";

    case "variant"
        if vmap.ok
            v_neg = map_variant(Tneg.wave_id, vmap);
            v_tgc = map_variant(Ttgc.wave_id, vmap);
            g_neg = v_neg; g_tgc = v_tgc;
            grp_labels = union_unique(v_neg, v_tgc);
        else
            warning('GroupBy=variant 但缺少 variant 映射，降级为 ALL');
            g_neg = repmat("ALL", height(Tneg), 1);
            g_tgc = repmat("ALL", height(Ttgc), 1);
            grp_labels = "ALL";
        end
        ytitle = "VARIANT";

    otherwise % 'variant_cluster'
        if vmap.ok
            v_neg = map_variant(Tneg.wave_id, vmap);
            v_tgc = map_variant(Ttgc.wave_id, vmap);
            c_neg = label_cluster(Tneg.cluster);
            c_tgc = label_cluster(Ttgc.cluster);
            % 组合
            g_neg = comb(vc(v_neg), vc(c_neg));
            g_tgc = comb(vc(v_tgc), vc(c_tgc));
            % 过滤 NA
            g_neg(g_neg=="NA") = [];
            g_tgc(g_tgc=="NA") = [];
            grp_labels = union_unique(g_neg, g_tgc);
            if isempty(grp_labels), grp_labels = "ALL"; g_neg(:)="ALL"; g_tgc(:)="ALL"; end
        else
            warning('GroupBy=variant_cluster 但缺少 variant 映射，降级为 cluster');
            g_neg = label_cluster(Tneg.cluster);
            g_tgc = label_cluster(Ttgc.cluster);
            grp_labels = union_unique(g_neg, g_tgc);
        end
        ytitle = "VARIANT | CLUSTER";
end
% 若标签为空，兜底为 ALL
if isempty(grp_labels)
    g_neg = repmat("ALL", height(Tneg), 1);
    g_tgc = repmat("ALL", height(Ttgc), 1);
    grp_labels = "ALL";
end
end

function s = vc(x), s = string(x); end
function g = comb(v, c)
ok = (v~="") & (c~="") & (c~="CNaN");
g = strings(numel(v),1); g(~ok) = "NA"; g(ok) = v(ok) + " | " + c(ok);
end
function lab = label_cluster(cnum)
% cnum 可能是 double（含 NaN）
cnum = double(cnum);
lab = "C" + string(cnum);
lab(~isfinite(cnum)) = "CNaN";
end
function u = union_unique(a,b)
u = unique([string(a(:)); string(b(:))], 'stable');
u(u=="NA" | u=="CNaN" | u=="") = [];
end
function v = map_variant(wave_id, vmap)
wid = string(wave_id);
[tf, loc] = ismember(wid, vmap.wave);
v = repmat("", numel(wid), 1);
v(tf) = vmap.var(loc(tf));
v(~tf) = "Unknown";
end

% ========== 组×参数 中位矩阵 ==========
function M = build_median_matrix_core(T, xname, xvec, g_for_rows, grp_labels)
xval = double(T.(xname));
y    = double(T.avoided_pct);
G = numel(grp_labels); nx = numel(xvec);
M = nan(G, nx);
[~,loc] = ismembertol(xval, xvec, 1e-12, 'DataScale', 1, 'ByRows', false);
for gi = 1:G
    gname = grp_labels(gi);
    idx_g = (g_for_rows == gname) & (loc>0);
    if ~any(idx_g), continue; end
    for j = 1:nx
        hit = idx_g & (loc == j);
        if any(hit), M(gi,j) = median(y(hit), 'omitnan'); end
    end
end
end

% ========== 数值标注（稀疏/限量） ==========
function local_annotate_heat(ax, M, xvec, yvec, fmt, fs, topK, qtile)
if nargin<5||isempty(fmt), fmt='%.1f'; end
if nargin<6||isempty(fs),  fs = 9;    end
if nargin<7||isempty(topK), topK = 0; end
if nargin<8||isempty(qtile), qtile = 1.0; end
xvec = double(xvec(:)).'; yvec = double(yvec(:)).';
[ny,nx] = size(M);
if numel(xvec)~=nx || numel(yvec)~=ny
    error('标注尺寸不匹配：numel(xvec)=%d, nx=%d;  numel(yvec)=%d, ny=%d', numel(xvec),nx,numel(yvec),ny);
end
MAX_CELLS_ANNOTATE = 3000;
numCells = nnz(isfinite(M));
if numCells > MAX_CELLS_ANNOTATE && topK==0 && qtile>=1.0
    qtile = 0.99; % 自动限流
end
mask = isfinite(M);
vals = abs(M(mask)); if isempty(vals), return; end
th = quantile(vals, qtile);
sel = false(size(M));
sel(mask & abs(M)>=th) = true;
if topK>0
    [~,I] = maxk(abs(M(sel)), min(topK, nnz(sel)));
    idx = find(sel);
    keep = false(size(sel));
    keep(idx(I)) = true;
    sel = keep;
end
axes(ax); hold(ax,'on');
cmap = colormap(ax); vlim = caxis(ax); nC = size(cmap,1);
[iinds, jinds] = find(sel);
for t=1:numel(iinds)
    i = iinds(t); j = jinds(t);
    v = M(i,j);
    tt = (v - vlim(1)) / max(eps, vlim(2)-vlim(1)); tt = min(max(tt,0),1);
    k = 1 + floor(tt*(nC-1)); rgb = cmap(k,:);
    L = 0.299*rgb(1) + 0.587*rgb(2) + 0.114*rgb(3);
    txtcol = [0 0 0]; if L < 0.55, txtcol = [1 1 1]; end
    text(xvec(j), yvec(i), sprintf(char(fmt), v), ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
        'FontSize', fs, 'FontWeight','bold', 'Color', txtcol, ...
        'Clipping','on','HitTest','off','PickableParts','none');
end
hold(ax,'off');
end

% ========== 学界常用高端色带 ==========
function cmap = get_nature_cmap(name, n)
if nargin<2||isempty(n), n = 256; end
name = lower(string(name));
switch name
    case "rdbu"   % 发散：蓝-白-红（零居中）
        c_lo = [0.230, 0.299, 0.754]; c_md = [1 1 1]; c_hi = [0.706, 0.016, 0.150];
        t = linspace(0,1,n).'; cmap = zeros(n,3);
        for i=1:n
            if t(i) < 0.5, a = t(i)/0.5;  cmap(i,:) = (1-a)*c_lo + a*c_md;
            else,          a = (t(i)-0.5)/0.5; cmap(i,:) = (1-a)*c_md + a*c_hi; end
        end
    case "viridis"
        base = [68 1 84; 71 44 122; 59 81 139; 44 113 142; 33 144 141; 39 173 129; 92 200 99; 170 220 50; 253 231 37]/255;
        cmap = interp1(linspace(0,1,size(base,1)), base, linspace(0,1,n));
    case "cividis"
        base = [0 32 76; 33 62 114; 60 89 124; 86 115 126; 112 142 128; 140 170 131; 172 199 133; 206 226 138; 236 250 143]/255;
        cmap = interp1(linspace(0,1,size(base,1)), base, linspace(0,1,n));
    case "inferno"
        base = [0 0 4; 31 12 72; 85 15 109; 136 34 106; 186 54 85; 227 89 51; 249 140 10; 252 195 65; 252 255 164]/255;
        cmap = interp1(linspace(0,1,size(base,1)), base, linspace(0,1,n));
    case "plasma"
        base = [13 8 135; 84 2 163; 139 10 165; 187 55 132; 225 101 97; 252 146 65; 254 195 40; 240 249 33]/255;
        cmap = interp1(linspace(0,1,size(base,1)), base, linspace(0,1,n));
    case "magma"
        base = [0 0 4; 28 16 68; 79 18 123; 129 37 129; 181 54 122; 229 80 100; 252 130 73; 253 183 98; 252 253 191]/255;
        cmap = interp1(linspace(0,1,size(base,1)), base, linspace(0,1,n));
    otherwise
        try, cmap = feval(char(name), n); catch, cmap = parula(n); end
end
end

% ========== 列名模糊检索 ==========
function out = detect_alt_col(T, keyword)
vnames = string(T.Properties.VariableNames);
idx = find(contains(vnames, string(keyword), 'IgnoreCase', true), 1);
if isempty(idx), out = ''; else, out = char(vnames(idx)); end
end


function labs = apply_cluster_order(labs, order, groupBy)
    labs = string(labs(:));
    gb = lower(string(groupBy));
    if isempty(labs), return; end

    % 显式顺序：ClusterOrder=[1 3 2 4] → "C1","C3","C2","C4"
    if isnumeric(order) && ~isempty(order)
        seq = "C" + string(order(:)');
        labs = seq(ismember(seq, labs));
        return
    end

    ord = lower(string(order));
    if ~ismember(ord,["ascend","asc","descend","desc"])
        return  % 默认保持 stable
    end

    if gb=="cluster"
        % 只含 "Ck"
        tok = regexp(labs,"^C(\d+)$",'tokens','once');
        num = cellfun(@(t) iff(isempty(t), NaN, str2double(t{1})), tok);
        [~,idx] = sort(num, ord);
        labs = labs(idx);

    elseif gb=="variant_cluster"
        % 形如 "Variant | Ck"：先按 Variant 字母序，再按 k
        var = regexprep(labs,"\s*\|\s*C\d+\s*$","");
        tok = regexp(labs,"C(\d+)$",'tokens','once');
        num = cellfun(@(t) iff(isempty(t), Inf, str2double(t{1})), tok); % 无数字放后
        T = table(string(var), num(:), (1:numel(labs))',...
                  'VariableNames',{'variant','cnum','idx'});
        T = sortrows(T, {'variant','cnum'}, {'ascend', ord});
        labs = labs(T.idx);
    end
end

function y = iff(cond, a, b), if cond, y=a; else, y=b; end, end
