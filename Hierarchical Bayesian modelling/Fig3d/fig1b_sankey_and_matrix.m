function fig = fig1b_sankey_and_matrix(results_root, varargin)
% Group→Target 桑基 + Signed bubble 矩阵（无缝节点、无描边气泡、紧凑图例）
% 口径：半标准化 + L1后验期望；颜色 S=(Pos-Neg)/W；气泡=同号后验份额
% 归一化：'share'(默认) 或 'none'
p = inputParser;
addParameter(p,'normalize_by_target','share',@(s) any(strcmpi(s,{'share','none'})));
addParameter(p,'use_posterior_expectation',true,@islogical);
addParameter(p,'standardize_to_sd',true,@islogical);
addParameter(p,'exclude_groups',["Region"],@(x)isstring(x)||iscellstr(x));
addParameter(p,'dpi',300,@(x)isnumeric(x)&&x>=72);
addParameter(p,'alpha_band',0.55,@(x)isnumeric(x));
addParameter(p,'curve_bias',[0.33 0.67],@(v)isnumeric(v)&&numel(v)==2);
addParameter(p,'palette_key','mainA8_deep',@(s)ischar(s)||isstring(s));
addParameter(p,'div_palette',[],@(x)isnumeric(x)&&size(x,2)==3||isempty(x));
addParameter(p,'matrix_size_exp',0.85,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'matrix_bubble_min_pts2',40,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'matrix_label_mode','auto',@(s)any(strcmpi(s,{'auto','all','none'})));
addParameter(p,'matrix_label_thresh',0.5,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);

parse(p,varargin{:}); A = p.Results;

if nargin<1||isempty(results_root), results_root='results_hbm'; end
outdir = fullfile(results_root,'M0_main','_figs'); if ~isfolder(outdir), mkdir(outdir); end

% ---- 分组映射（默认排除 Region）----
group_map = struct('Population',["age" "student" "education"], ...
    'Healthcare',["hospital" "beds"], 'Economy',["gdp" "income"], ...
    'Urbanicity',["urban" "density" "tertiary"], 'Mobility',["passenger" "transit"], ...
    'Tourism',["inbound" "overnight" "visitors"], 'Region',["region"], 'Climate',["climate"]);
for k=1:numel(A.exclude_groups)
    if isfield(group_map, char(A.exclude_groups(k))), group_map = rmfield(group_map, char(A.exclude_groups(k))); end
end

% ---- 载入缓存 ----
C = load(fullfile(results_root,'cache','hbm_cache.mat')); cache=C.cache;
targets = string(fieldnames(cache.Y)); fe_names = string(cache.fe_names_orig);
Xfe = cache.X_fe; Kc = numel(cache.cont_orig);

% 目标 SD(Y)
sdY = zeros(numel(targets),1);
for ti=1:numel(targets)
    yv = cache.Y.(targets(ti)); ok = isfinite(yv) & all(isfinite(Xfe),2);
    sdy = std(yv(ok),'omitnan'); if ~isfinite(sdy)||sdy<=0, sdy=1; end, sdY(ti)=sdy;
end

% ---- 组×目标矩阵 ----
gname = string(fieldnames(group_map)); G=numel(gname); T=numel(targets);
Pos=zeros(G,T); Neg=zeros(G,T); Sig=zeros(G,T);
phi=@(z)exp(-0.5*z.^2)/sqrt(2*pi);

for ti=1:T
    Tab = readtable(fullfile(results_root,'M0_main',targets(ti),'summary.csv'),'TextType','string');
    is_named = startsWith(Tab.param,"beta:"); is_idx = ~is_named & startsWith(Tab.param,"beta_");
    A1=table(); B1=table();
    if any(is_named)
        tmp=Tab(is_named,:); A1=table(extractAfter(tmp.param,"beta:"),tmp.mean,tmp.sd, ...
            'VariableNames',{'param','mean','sd'});
    end
    if any(is_idx)
        tmp=Tab(is_idx,:); k=double(extractAfter(tmp.param,"beta_"));
        nm=strings(height(tmp),1);
        for r=1:height(tmp)
            if ~isnan(k(r)) && k(r)>=1 && k(r)<=numel(fe_names), nm(r)=fe_names(k(r));
            else, nm(r)=string(tmp.param(r)); end
        end
        B1=table(nm,tmp.mean,tmp.sd,'VariableNames',{'param','mean','sd'});
    end
    Tb=[A1;B1]; if isempty(Tb), continue; end

    % 半标准化
    idx_in_fe = nan(height(Tb),1);
    for r=1:height(Tb), hit=find(fe_names==Tb.param(r),1,'first'); if ~isempty(hit), idx_in_fe(r)=hit; end, end
    sdx = ones(height(Tb),1);
    okrows = isfinite(cache.Y.(targets(ti))) & all(isfinite(Xfe),2);
    for r=1:height(Tb)
        j=idx_in_fe(r);
        if ~isnan(j)
            if j<=Kc, sdx(r)=1;
            else, p1=mean(Xfe(okrows,j)>0.5,'omitnan'); sdx(r)=sqrt(max(p1*(1-p1),eps));
            end
        end
    end
    mu0=Tb.mean; sd0=Tb.sd; sd0(~isfinite(sd0)|sd0<=0)=eps;
    if A.standardize_to_sd, mu=mu0.*(sdx./max(sdY(ti),eps)); sd=sd0.*(sdx./max(sdY(ti),eps));
    else, mu=mu0; sd=sd0; end

    for g=1:G
        toks=string(group_map.(gname(g)));
        m=false(height(Tb),1); for tt=1:numel(toks), m=m|contains(lower(Tb.param),lower(toks(tt))); end
        if ~any(m), continue; end
        mu_g=mu(m); sd_g=sd(m);
        a=mu_g./sd_g;
        Epos = sd_g.*phi(a) + mu_g.*normcdf(a,0,1);
        Eneg = sd_g.*phi(a) - mu_g.*(1-normcdf(a,0,1));
        Pos(g,ti)=nansum(max(Epos,0)); Neg(g,ti)=nansum(max(Eneg,0));
        z=abs(mu_g./sd_g); p_same=normcdf(z,0,1); Sig(g,ti)=mean(p_same>=0.95,'omitnan');
    end
end

W=Pos+Neg; S=(Pos-Neg)./max(W,eps);



% ======== NEW: export domain-level strengths & direction balance ========
% 保留一份“原始顺序”的 W/S 以便导出（不受后续排序影响）
W_raw = W;                % G×T
S_raw = S;                % G×T
Sig_raw = Sig;            % G×T

% 计算“按目标归一化”的份额（share-by-target）
colsum = sum(W_raw, 1);                     % 1×T
colsum(colsum<=0) = eps;
Share_by_target = W_raw ./ colsum(ones(size(W_raw,1),1), :);  % G×T，矩阵广播

% 组名、目标名
Gnames = string(gname(:));                  % G×1
Tnames = string(targets(:));                % T×1

% ------- 组织长表（Target, Group, W, S, Share, Sig） -------
G = numel(Gnames); T = numel(Tnames);
ng = G*T;
Target_col = strings(ng,1);
Group_col  = strings(ng,1);
W_col      = zeros(ng,1);
S_col      = zeros(ng,1);
Share_col  = zeros(ng,1);
Sig_col    = zeros(ng,1);

% 用矩阵展开到长向量（不使用动态点索引）
% 说明：按列主序展开 -> 先沿组方向，再推进目标
W_vec = reshape(W_raw,   G*T, 1);
S_vec = reshape(S_raw,   G*T, 1);
Q_vec = reshape(Share_by_target, G*T, 1);
P_vec = reshape(Sig_raw, G*T, 1);

% 目标与组的配对（同样矩阵式赋值）
% 生成每个目标重复 G 次的向量、每个组重复 T 次的向量
Target_rep = repmat(Tnames', G, 1);   % G×T
Group_rep  = repmat(Gnames,  1, T);   % G×T

Target_col(:) = reshape(Target_rep, G*T, 1);
Group_col(:)  = reshape(Group_rep,  G*T, 1);
W_col(:)      = W_vec;
S_col(:)      = S_vec;
Share_col(:)  = Q_vec;
Sig_col(:)    = P_vec;

T_long = table(cellstr(Target_col), cellstr(Group_col), ...
               W_col, S_col, Share_col, Sig_col, ...
               'VariableNames', {'Target','Group','W','S','ShareWithinTarget','SigShare'});

% 保存 CSV
outdir = fullfile(results_root,'M0_main','_figs');
if ~isfolder(outdir), mkdir(outdir); end
writetable(T_long, fullfile(outdir,'Fig1b_domain_strengths_long.csv'));

% 额外：保存“宽表”版本（每个表一类量，行=Group，列=Target）
% 为了合法列名，转成有效变量名（不改变文件内显示）
colnames = matlab.lang.makeValidName(cellstr(Tnames));
T_W_wide = array2table(W_raw,   'VariableNames', colnames, 'RowNames', cellstr(Gnames));
T_S_wide = array2table(S_raw,   'VariableNames', colnames, 'RowNames', cellstr(Gnames));
T_Q_wide = array2table(Share_by_target, 'VariableNames', colnames, 'RowNames', cellstr(Gnames));

writetable(T_W_wide, fullfile(outdir,'Fig1b_W_wide.csv'), 'WriteRowNames',true);
writetable(T_S_wide, fullfile(outdir,'Fig1b_S_wide.csv'), 'WriteRowNames',true);
writetable(T_Q_wide, fullfile(outdir,'Fig1b_W_share_by_target_wide.csv'), 'WriteRowNames',true);

% 命令行打印：每个目标按 W 的份额（share-by-target）降序显示前 3 个组
fprintf('[Fig1b] Domain strength W and direction S exported.\n');
for t = 1:T
    [~, idx] = sort(Share_by_target(:,t), 'descend');
    topk = idx(1:min(3,G));
    fprintf('  Target %-12s | ', char(Tnames(t)));
    for k = 1:numel(topk)
        g = topk(k);
        fprintf('%s: share=%.2f  S=%.2f  ', char(Gnames(g)), Share_by_target(g,t), S_raw(g,t));
    end
    fprintf('\n');
end
% ======== NEW block end ========






switch lower(A.normalize_by_target)
    case 'share', colsum=sum(W,1); colsum(colsum==0)=1; W_plot=W./colsum; subtitle='share within target';
    otherwise,    W_plot=W; subtitle='absolute (semi-standardized)';
end
absSum=sum(W,1);

% ---- 配色 ----
if isempty(A.div_palette), negC=[0.10 0.35 0.70]; posC=[0.82 0.10 0.12];
else, negC=A.div_palette(1,:); posC=A.div_palette(end,:); end
neuC=[0.82 0.82 0.82];
blendDiv=@(s)(s>=0).*(neuC+(posC-neuC).*abs(s))+(s<0).*(neuC+(negC-neuC).*abs(s));

try, Cg = smf_colorlib(A.palette_key); catch, Cg = smf_colorlib('classic10'); end
if size(Cg,1)<G, Cg=repmat(Cg,ceil(G/size(Cg,1)),1); end, Cg=Cg(1:G,:);

% ---- 画布 ----
fig=figure('Color','w','Units','pixels','Position',[60 60 1550 920]);
axL=axes('Parent',fig,'Units','normalized','Position',[0.07 0.15 0.60 0.77]); hold(axL,'on');
set(axL,'XLim',[0 1],'YLim',[0 1],'Box','off','Layer','top','XColor','none');
axR=axes('Parent',fig,'Units','normalized','Position',[0.72 0.15 0.24 0.77]); hold(axR,'on');

% === 几何（无缝对齐）===
y_margin=0.04; flow_ylo=y_margin; flow_yhi=1-y_margin; gap_vert=0.015;
node_w=0.060;                         % 节点宽
xL_node=0.08;                         % 左节点 x
xL_flow0=xL_node+node_w;              % 流起点紧贴左节点右边缘
xR_flow1=0.86;                        % 流终点
xR_node=xR_flow1;                     % 右节点与流终点紧贴

[left_blocks,grp_centers]=stack_blocks_range(sum(W_plot,2),gap_vert,flow_ylo,flow_yhi);
[right_blocks,tgt_centers]=stack_blocks_range(sum(W_plot,1)',gap_vert,flow_ylo,flow_yhi);
[Lseg,Rseg,CL,CR]=stack_links(W_plot,left_blocks,right_blocks);

% 流（由细到粗）
edge_vals=W_plot; edge_vals(edge_vals<=0)=nan;
[~,ord]=sort(edge_vals(:),'ascend');
nbezier=28; c1=A.curve_bias(1); c2=A.curve_bias(2);
for k=1:numel(ord)
    id=ord(k); if isnan(id), continue; end
    [g,t]=ind2sub([G T],id); if W_plot(g,t)<=0, continue; end
    yL0=Lseg(g).y0(t); yL1=Lseg(g).y1(t);
    yR0=Rseg(t).y0(g); yR1=Rseg(t).y1(g);
    [xb_up,yb_up]=bezier_band([xL_flow0,yL1],[xR_flow1,yR1],nbezier,c1,c2);
    [xb_lo,yb_lo]=bezier_band([xR_flow1,yR0],[xL_flow0,yL0],nbezier,c1,c2);
    col=blendDiv(S(g,t));
    h=patch(axL,'XData',[xb_up(:);xb_lo(:)],'YData',[yb_up(:);yb_lo(:)], ...
        'FaceColor',col,'FaceAlpha',A.alpha_band,'EdgeColor','none'); safe_unpick(h);
end

% 左侧组名刻度
set(axL,'YTick',grp_centers,'YTickLabel',gname,'YColor',[0.25 0.25 0.25], ...
    'FontSize',12,'TickDir','out','TickLength',[0.004 0.004]);
line(axL,[xL_flow0 xL_flow0],[flow_ylo flow_yhi],'Color',[0.88 0.88 0.88],'LineWidth',0.6,'HitTest','off');

% 左右节点（无缝贴合）
for g=1:G
    b=left_blocks(g);
    h=roundedRect(axL,xL_node,b.y0,node_w,b.y1-b.y0,0.006,Cg(g,:),[0 0 0]*0.55,0.8,0.98); safe_unpick(h);
end
for t=1:T
    b=right_blocks(t);
    h=roundedRect(axL,xR_node-node_w,b.y0,node_w,b.y1-b.y0,0.006,[0.15 0.15 0.15],[0 0 0]*0.6,0.8,1.0); safe_unpick(h);
    ht=text(axL,xR_node+0.012,tgt_centers(t),sprintf('%s',targets(t)), ...
        'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',12,'Interpreter','none'); safe_pick(ht);
    ht2=text(axL,xR_node+0.012,tgt_centers(t)-0.03,sprintf('\\Sigma=%.2f',absSum(t)), ...
        'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',9,'Color',[0.35 0.35 0.35],'Interpreter','tex'); safe_pick(ht2);
end

% 中线气泡（无描边）
Yc_raw=0.5*(CL+CR); Yc_nudged=Yc_raw; minGap=0.012;
for t=1:T, Yc_nudged(:,t)=nudge_y(Yc_raw(:,t),minGap); end
for t=1:T
    for g=1:G
        if W_plot(g,t)<=0, continue; end
        sig=max(0,min(1,Sig(g,t)));
        s=70+(900-70)*(sig^1.6);
        jitter=(( -1)^g)*0.010*(0.5+0.5*rand());
        xc=0.5*(xL_flow0+xR_flow1)+jitter; yc=Yc_nudged(g,t);
        bcol=blendDiv(S(g,t));
        h=scatter(axL,xc,yc,s,'o','filled','MarkerFaceColor',bcol, ...
            'MarkerEdgeColor','none','MarkerFaceAlpha',0.95); safe_pick(h);
        show_label = false;
        switch lower(A.matrix_label_mode)
            case 'all'
                show_label = true;
            case 'none'
                show_label = false;
            otherwise % 'auto'
                show_label = (Sig(g,t) >= A.matrix_label_thresh);
        end
        if show_label
            ht = text(axR, t, g, sprintf('%d%%', round(100*Sig(g,t))), ...
                'HorizontalAlignment','center','VerticalAlignment','middle', ...
                'FontSize',9,'Color',[0.05 0.05 0.05]);
            safe_pick(ht);
        end

    end
end
safe_pick(title(axL,sprintf('Group \\rightarrow Target Sankey with significance bubbles (%s)',subtitle), ...
    'FontSize',18,'FontWeight','bold','Interpreter','none'));

% === 紧凑可选中的气泡图例（独立小轴，可在编辑器中拖动） ===
leg_pos=[0.34 0.06 0.32 0.07];
legend_ax = axes('Parent',fig,'Position',leg_pos,'Color','none','XTick',[],'YTick',[]);
hold(legend_ax,'on'); axis(legend_ax,'off');
set(legend_ax,'PickableParts','all','HitTest','on');   % 可选中移动
rep=[0 0.33 0.66 1.0]; labs={'0%','33%','66%','100%'}; xs=linspace(0.10,0.90,numel(rep));
for i=1:numel(rep)
    si = 70+(900-70)*(rep(i)^1.6);
    h=scatter(legend_ax,xs(i),0.5,si,'filled','MarkerFaceColor',[0.7 0.7 0.7], ...
        'MarkerEdgeColor','none','MarkerFaceAlpha',0.95); safe_pick(h);
    ht=text(legend_ax,xs(i),0.05,labs{i},'HorizontalAlignment','center','FontSize',10); safe_pick(ht);
end
% ht=text(legend_ax,0.0,0.92,'Bubble = share with posterior same-sign prob (P_{same}) \ge 0.95', ...
%     'HorizontalAlignment','left','FontSize',10,'Interpreter','tex'); safe_pick(ht);

text(legend_ax, 0.0, 0.92, ...
    'Bubble = share with posterior same-sign prob (P_{same} \geq 0.95)', ...
    'HorizontalAlignment','left','FontSize',10,'Interpreter','tex');


% === 右：Signed bubble matrix（无描边 + 自适应） ===
imagesc(axR,[1 T],[1 G], zeros(G,T)); axis(axR,'ij'); colormap(axR, flipud(gray(2))); caxis(axR,[0 1]);
set(axR,'XTick',1:T,'XTickLabel',targets,'XTickLabelRotation',45,'YTick',1:G,'YTickLabel',gname,'FontSize',11,'TickDir','out');
grid(axR,'on'); axR.GridAlpha=0.12;

drawnow; pp_ax=getpixelposition(axR); ppi=get(groot,'ScreenPixelsPerInch');
cell_px=min(pp_ax(3)/max(T,1.2), pp_ax(4)/max(G,1.2)); cell_pts=cell_px*72/ppi;
smax=(0.58*cell_pts)^2; smin=A.matrix_bubble_min_pts2; wmax=max(W_plot(:)); if wmax<=0, wmax=1; end
for g=1:G
    for t=1:T
        w=W_plot(g,t); if w<=0, continue; end
        col=blendDiv(S(g,t));
        sz=smin+(smax-smin)*(w/wmax)^A.matrix_size_exp;
        h=scatter(axR,t,g,sz,'o','filled','MarkerFaceColor',col,'MarkerEdgeColor','none','MarkerFaceAlpha',0.95); safe_pick(h);
        if Sig(g,t)>=0.5
            ht=text(axR,t,g,sprintf('%d%%',round(100*Sig(g,t))), ...
                'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',9,'Color',[0.05 0.05 0.05]); safe_pick(ht);
        end
    end
end
safe_pick(title(axR,'Signed bubble matrix (same scale & color as Sankey)','FontSize',14,'FontWeight','bold','Interpreter','none'));

% 色条（S）
cbax=axes('Position',[0.965 0.15 0.015 0.77]); axis(cbax,'off');
Ncb=64; Scb=linspace(-1,1,Ncb)'; Ccb=zeros(Ncb,3); for i=1:Ncb, Ccb(i,:)=blendDiv(Scb(i)); end
image(cbax,'XData',[0 1],'YData',[0 1],'CData',permute(Ccb,[1 3 2])); set(cbax,'YDir','normal');
text(cbax,1.4,1,'S=+1','FontSize',9,'HorizontalAlignment','left');
text(cbax,1.4,0.5,'0','FontSize',9,'HorizontalAlignment','left');
text(cbax,1.4,0,'S=-1','FontSize',9,'HorizontalAlignment','left');

% 顶部注释
annotation(fig,'textbox',[0.06 0.935 0.88 0.05], ...
    'String','Fig. 1b | Group-to-target influence (Sankey + signed bubble matrix)', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'EdgeColor','none','FontSize',16,'FontWeight','bold','Interpreter','none');

boxed_note(fig,[0.065 0.91 0.40 0.06], ...
    "Semi-standardized L1 strength" + newline + ...
    "\betã = \beta·SD(X)/SD(Y_{target});  one-hot SD(X)=\surd(p(1-p))" + newline + ...
    "Width = Pos+Neg;  Color S=(Pos-Neg)/W;  " + subtitle);

make_editable(fig);
print(fig, fullfile(outdir, sprintf('Fig1b_sankey_matrix_%s.png',lower(A.normalize_by_target))), ...
    '-dpng', sprintf('-r%d',A.dpi));
end

% ===== 工具 =====
function [blocks, centers]=stack_blocks_range(weights,gap,ylo,yhi)
w=max(weights(:),0); S=sum(w); if S<=0, w(:)=1; S=sum(w); end
G=numel(w); usable=max(yhi-ylo-gap*(G-1),1e-6);
h=usable*(w/S); y=zeros(G,2); y0=ylo;
for i=1:G, y1=y0+h(i); y(i,:)=[y0 y1]; y0=y1+gap; end
blocks=struct('y0',num2cell(y(:,1)),'y1',num2cell(y(:,2))); centers=mean(y,2);
end
function [Lseg,Rseg,CL,CR]=stack_links(Wn,Lblocks,Rblocks)
[G,T]=size(Wn);
Lseg(G).y0=[]; Lseg(G).y1=[]; Rseg(T).y0=[]; Rseg(T).y1=[];
for g=1:G
    b=Lblocks(g); Hg=b.y1-b.y0; share=Wn(g,:); s=sum(share); if s<=0, share(:)=0; s=1; end
    hh=Hg*(share/s); y0=b.y0; y0s=zeros(1,T); y1s=zeros(1,T);
    for t=1:T, y1=y0+hh(t); y0s(t)=y0; y1s(t)=y1; CL(g,t)=(y0+y1)/2; y0=y1; end
    Lseg(g).y0=y0s; Lseg(g).y1=y1s;
end
for t=1:T
    b=Rblocks(t); Ht=b.y1-b.y0; share=Wn(:,t); s=sum(share); if s<=0, share(:)=0; s=1; end
    hh=Ht*(share/s); y0=b.y0; y0s=zeros(1,G); y1s=zeros(1,G);
    for g=1:G, y1=y0+hh(g); y0s(g)=y0; y1s(g)=y1; CR(g,t)=(y0+y1)/2; y0=y1; end
    Rseg(t).y0=y0s; Rseg(t).y1=y1s;
end
end
function [xb,yb]=bezier_band(p0,p1,n,c1,c2)
x0=p0(1); y0=p0(2); x1=p1(1); y1=p1(2); dx=x1-x0;
C1=[x0+c1*dx, y0]; C2=[x0+c2*dx, y1];
t=linspace(0,1,n)'; B0=(1-t).^3; B1=3*(1-t).^2.*t; B2=3*(1-t).*t.^2; B3=t.^3;
xb=B0*x0 + B1*C1(1) + B2*C2(1) + B3*x1;
yb=B0*y0 + B1*C1(2) + B2*C2(2) + B3*y1;
end
function y=nudge_y(y, minGap)
[ys,ord]=sort(y);
for i=2:numel(ys), if ys(i)-ys(i-1)<minGap, ys(i)=ys(i-1)+minGap; end, end
ys=max(1e-4,min(1-1e-4,ys)); y(ord)=ys;
end

function h = roundedRect(ax, x, y, w, hgt, r, faceColor, edgeColor, lw, alpha)
% 圆角矩形（无缝贴合流带）— 使用 FILL，避免 patch 的第4参数尺寸限制
r = max(0, min([r, w/2, hgt/2]));
n = 24;
th1 = linspace(pi,   1.5*pi, n);   % 左下
th2 = linspace(1.5*pi, 2*pi, n);   % 右下
th3 = linspace(0,   0.5*pi, n);    % 右上
th4 = linspace(0.5*pi,   pi, n);   % 左上

xb = [ (x+r)     + r*cos(th1), ...
    (x+w-r)   + r*cos(th2), ...
    (x+w-r)   + r*cos(th3), ...
    (x+r)     + r*cos(th4) ];

yb = [ (y+r)     + r*sin(th1), ...
    (y+r)     + r*sin(th2), ...
    (y+hgt-r) + r*sin(th3), ...
    (y+hgt-r) + r*sin(th4) ];

h = fill(ax, xb, yb, faceColor, ...
    'EdgeColor', edgeColor, 'LineWidth', lw, 'FaceAlpha', alpha, ...
    'HitTest','off','PickableParts','none');
end




function boxed_note(fig,pos,str)
hta=annotation(fig,'textbox',pos,'String',str,'EdgeColor',[0.85 0.85 0.85], ...
    'BackgroundColor',[1 1 1],'FitBoxToText','off','FontSize',10, ...
    'HorizontalAlignment','left','VerticalAlignment','middle','Interpreter','tex');
set(hta,'PickableParts','all','HitTest','on');
end
function make_editable(fig)
set(findall(fig,'Type','axes'),'Layer','top');
set(findall(fig,'Type','text'),'HitTest','on');
set(findall(fig,'Type','scatter'),'HitTest','on');
hP=findall(fig,'Type','patch'); for k=1:numel(hP), if isprop(hP(k),'HitTest'), set(hP(k),'HitTest','off'); end, end
end

function safe_unpick(h)
% 让对象在图形编辑器里“不可选中”（不抢点击）
try
    if isgraphics(h)
        if isprop(h,'HitTest'),        set(h,'HitTest','off');        end
        if isprop(h,'PickableParts'),  set(h,'PickableParts','none');  end
    end
catch
end
end

function safe_pick(h)
% 让对象在图形编辑器里可选中（标题、标签、散点、图例文本等）
try
    if isgraphics(h)
        if isprop(h,'HitTest'),        set(h,'HitTest','on');         end
        if isprop(h,'PickableParts'),  set(h,'PickableParts','all');   end
    end
catch
end
end
