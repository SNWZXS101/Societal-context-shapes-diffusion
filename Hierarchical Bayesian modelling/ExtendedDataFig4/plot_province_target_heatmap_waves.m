function H = plot_province_target_heatmap_waves(varargin)
% 简版：30省 × 目标 的省级随机效应热图（WAVES 权重，逐格patch，可编辑）
% 依赖：smf_colorlib.m；需要 results_hbm/* 与 features_aug.mat
%
% 例：
% H = plot_province_target_heatmap_waves_simple( ...
%       'results_root','results_hbm', ...
%       'features_mat','features_aug.mat', ...
%       'primaryTarget','delta_early', ...   % 若 sortBy='primaryDesc' 时用于排序
%       'sortBy','none', ...                 % 'none' | 'primaryDesc' | 'alpha'
%       'FixedClim',[-0.2 0.2], ...          % 固定色域；留空则用对称 p98
%       'cmap','icefire', ...                % 来自 smf_colorlib
%       'labelFontSize',9, ...
%       'title','Province × target heatmap (WAVES weighting)');

%% ---- 参数 ----
p = inputParser;
p.addParameter('results_root','results_hbm',@(s)ischar(s)||isstring(s));
p.addParameter('features_mat','features_aug.mat',@(s)ischar(s)||isstring(s));
p.addParameter('primaryTarget','delta_early',@(s)ischar(s)||isstring(s));
p.addParameter('sortBy','none',@(s)any(strcmpi(s,{'none','primaryDesc','alpha'})));
p.addParameter('FixedClim',[],@(x)isnumeric(x)&&numel(x)==2 || isempty(x));
p.addParameter('cmap','icefire',@(s)ischar(s)||isstring(s));
p.addParameter('renderer','painters',@(s)ischar(s)||isstring(s));
p.addParameter('labelFontSize',9,@(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('title','Province × target heatmap (WAVES weighting)',@(s)ischar(s)||isstring(s));
p.parse(varargin{:});
args = p.Results;

fontMain = pick_font({'Arial','Helvetica','Nimbus Sans L'});
Cgray    = [0.2 0.2 0.2];

%% ---- 读 cache / 省 / 目标 ----
cache = load(fullfile(args.results_root,'cache','hbm_cache.mat')); cache = cache.cache;
prov_levels = string(cache.prov_names(:));
targets_all = string(cache.active_targets(:)');

%% ---- 读 features 并算 WAVES 权重（按波数归一） ----
Tfeat = load_features(args.features_mat);
Tfeat.variant = canon_variant(pick_col_val(Tfeat,["variant","Variant"]));
Tfeat.variant = replace(Tfeat.variant,"Original","Alpha"); % 与 ranef 对齐
var_levels_all = unique(Tfeat.variant,'stable');
W = waves_weights_by_count(Tfeat, prov_levels, var_levels_all);   % P×K

%% ---- 逐目标装配 u_i^adj ----
P = numel(prov_levels); Tn = numel(targets_all);
Uadj = nan(P, Tn);
for t = 1:Tn
    yname = targets_all(t);
    rp = readtable(fullfile(args.results_root,'M0_main',char(yname),'ranef_province.csv'));
    rv = readtable(fullfile(args.results_root,'M0_main',char(yname),'ranef_variant.csv'));
    rv.variant = canon_variant(string(rv.variant));
    rv.variant = replace(rv.variant,"Original","Alpha");

    mu_var = nan(numel(var_levels_all),1);
    for k=1:numel(var_levels_all)
        j = find(rv.variant==var_levels_all(k),1);
        if ~isempty(j), mu_var(k) = rv.mean(j); end
    end
    mu_prov = nan(P,1);
    for i=1:P
        j = find(string(rp.province)==prov_levels(i),1);
        if ~isempty(j), mu_prov(i) = rp.mean(j); end
    end
    Uadj(:,t) = mu_prov + W * mu_var;
end

%% ---- 排序（可选） ----
order_idx = (1:P)';
switch lower(args.sortBy)
    case 'primarydesc'
        [okPri, colPri] = ismember(lower(string(args.primaryTarget)), lower(targets_all));
        if ~okPri, colPri = 1; end
        [~, order_idx] = sort(Uadj(:,colPri), 'descend');
    case 'alpha'
        [~, order_idx] = sort(prov_levels);
    otherwise
        % 'none'：保持原始省顺序
end
Uadj_s = Uadj(order_idx,:); prov_s = prov_levels(order_idx);

%% ---- 颜色：发散，0居中；对称 p98 或固定 ----
vals = Uadj_s(~isnan(Uadj_s)); if isempty(vals), vals=0; end
if isempty(args.FixedClim)
    q = prctile(abs(vals),98); q = max(q, eps);
    clim = [-q q];
else
    clim = args.FixedClim(:).';
end
Cmap = smf_colorlib(args.cmap,'K',201);
mapfun = @(x) value2rgb(x, clim, Cmap);

%% ---- 手工布局：主热图 + 底部色条（逐格patch，可编辑） ----
fig = figure('Color','w','Name','Heatmap (WAVES)','Renderer',char(args.renderer));
posMain = [0.12 0.15 0.78 0.75];
posCB   = [0.12 0.09 0.78 0.04];

% 标题：annotation，可拖拽
annotation(fig,'textbox',[0.12 0.93 0.76 0.05], 'String',char(args.title), ...
    'EdgeColor','none','HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontName',fontMain,'FontSize',11,'FontWeight','normal','Color',Cgray, ...
    'Interpreter','none','HitTest','on','PickableParts','all');

% 主热图
ax = axes('Parent',fig,'Position',posMain,'Units','normalized'); hold(ax,'on');
[RR,CC] = size(Uadj_s);
for r=1:RR
    for c=1:CC
        v = Uadj_s(r,c);
        rgb = mapfun(v);
        x0=c-0.5; x1=c+0.5; y0=r-0.5; y1=r+0.5;
        patch('Parent',ax,'XData',[x0 x1 x1 x0],'YData',[y0 y0 y1 y1], ...
              'FaceColor',rgb,'EdgeColor',[0.94 0.94 0.94],'LineWidth',0.5, ...
              'HitTest','on','PickableParts','all');
    end
end
set(ax,'XLim',[0.5 CC+0.5],'YLim',[0.5 RR+0.5],'YDir','reverse', ...
    'XTick',1:CC,'XTickLabel',targets_all, ...
    'YTick',1:RR,'YTickLabel',prov_s, ...
    'TickDir','out','TickLength',[0.002 0.002], ...
    'FontName',fontMain,'FontSize',args.labelFontSize, ...
    'XColor',Cgray,'YColor',Cgray,'LineWidth',0.6,'Box','on');
xtickangle(ax,30); hold(ax,'off');

% —— 竖直色条（右侧，带数值刻度）——
posCB = [0.91 0.15 0.02 0.75];              % 右侧竖条位置
axCB = axes('Parent',fig,'Position',posCB,'Units','normalized');

hold(axCB,'on');  % 不要 axis off
K = size(Cmap,1);
for k = 1:K
    y0 = (k-1)/K; y1 = k/K; x0 = 0; x1 = 1;
    patch('Parent',axCB,'XData',[x0 x1 x1 x0],'YData',[y0 y0 y1 y1], ...
          'FaceColor',Cmap(k,:),'EdgeColor','none', ...
          'HitTest','on','PickableParts','all');
end

% 刻度与标签（数值来自 clim）
yticks = [0 0.25 0.5 0.75 1];
ytlabs  = arrayfun(@(z)sprintf('%.2g',z), linspace(clim(1),clim(2),numel(yticks)), 'uni',0);

set(axCB,'XLim',[0 1],'YLim',[0 1], ...
    'YTick',yticks,'YTickLabel',ytlabs, ...
    'XTick',[], ...
    'YAxisLocation','right', ...           % 刻度在右侧
    'TickDir','out','TickLength',[0.01 0.01], ...
    'FontName',fontMain,'FontSize',args.labelFontSize, ...
    'YColor',Cgray,'LineWidth',0.6, ...
    'Color','none','Box','on','Layer','top');   % 半透明背景、上层绘制

% 竖直标题（可选中）
text(axCB, 1.8, 0.5, 'Adjusted province RE (WAVES; 0-centered)', ...
    'Rotation',90,'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontName',fontMain,'FontSize',args.labelFontSize,'Color',Cgray, ...
    'HitTest','on','PickableParts','all');

hold(axCB,'off');



% 让图内元素可选
set(findall(fig,'-property','PickableParts'),'PickableParts','all');
set(findall(fig,'-property','HitTest'),'HitTest','on');

% 返回
H = struct('fig',fig,'ax',ax,'axCB',axCB,'order_idx',order_idx, ...
           'targets',targets_all,'clim',clim);

%% ===== helpers =====
function T = load_features(fpath)
    if exist(fpath,'file')~=2, error('未找到 features_mat：%s', fpath); end
    info = whos('-file', fpath);
    cand = ["Features","features","T","Table","tbl","Features_aug","features_aug"];
    names = string({info.name}); hit = "";
    for ii=1:numel(cand)
        j = find(strcmpi(names, cand(ii)), 1);
        if ~isempty(j), hit = names(j); break; end
    end
    if strlength(hit)>0
        S = load(fpath, char(hit)); T = S.(char(hit)); if istable(T), return; end
    end
    isTab = strcmp({info.class}, 'table');
    if any(isTab), nm = info(find(isTab,1)).name; S=load(fpath,nm); T=S.(nm); return; end
    error('Features 中未找到 table 变量：%s', fpath);
end

function v = pick_col_val(T, cands)
    vn = string(T.Properties.VariableNames);
    for i=1:numel(cands)
        j = find(strcmpi(vn, cands(i)), 1);
        if ~isempty(j), v = string(T.(vn(j))); return; end
    end
    error('Features 缺少列：%s', strjoin(cellstr(cands),', '));
end

function v = canon_variant(v)
    v = string(v); L = lower(strtrim(v)); v = strings(size(L));
    for i=1:numel(L)
        s=L(i);
        if any(strcmp(s,["original","ancestral","wuhan","wildtype","wild-type","wild type"])), v(i)="Alpha";
        elseif s=="alpha", v(i)="Alpha";
        elseif s=="delta", v(i)="Delta";
        elseif s=="omicron", v(i)="Omicron";
        else, v(i)=upper(extractBefore(L(i),1))+extractAfter(L(i),1);
        end
    end
end

function W = waves_weights_by_count(T, prov_levels, var_levels)
    prov_ser = pick_col_val(T,["province","Province"]);
    var_ser  = canon_variant(pick_col_val(T,["variant","Variant"]));
    var_ser  = replace(var_ser,"Original","Alpha");
    P = numel(prov_levels); K = numel(var_levels);
    W = zeros(P,K);
    for i=1:P
        mP = prov_ser==prov_levels(i);
        for k=1:K
            mV = var_ser==var_levels(k);
            W(i,k) = sum(mP & mV);  % 波数
        end
        s = sum(W(i,:)); if s>0, W(i,:)=W(i,:)/s; end
    end
end

function rgb = value2rgb(x,clim,Cmap)
    if ~isfinite(x), rgb=[0.95 0.95 0.95]; return; end
    lo=clim(1); hi=clim(2); if hi<=lo, hi=lo+eps; end
    t=(x-lo)/(hi-lo); t=min(max(t,0),1);
    idx=1+floor(t*(size(Cmap,1)-1)); rgb=Cmap(idx,:);
end

function fn = pick_font(cands)
    try
        avail=listfonts;
        for ii=1:numel(cands), if any(strcmpi(avail,cands{ii})), fn=cands{ii}; return; end, end
        fn=avail{1};
    catch, fn='Helvetica';
    end
end
end
