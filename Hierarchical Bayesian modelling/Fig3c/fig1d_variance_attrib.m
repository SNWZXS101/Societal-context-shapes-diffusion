function fig = fig1d_variance_attrib(results_root, varargin)
% FIG1D_VARIANCE_ATTRIB — 目标维度的方差归因（水平堆叠条 + 数值导出）
% 分量：FE / Province RE (含斜率) / Variant RE (含斜率) / Residual
% 尺度：变换后尺度。支持异方差与 Student-t 噪声。
%
% 额外输出（新增）：
%   results_root/M0_main/_figs/Fig1d_variance_shares.csv
%   每个目标：四分量的 posterior mean 与 95%CrI（2.5%,97.5%），并同时打印到命令行。
%
% 计算口径（每抽样 s，先算分量方差，再“先归一化后平均”）：
%   V_FE(s)   = Var_n( X*beta_s )
%   V_Prov(s) = Var_n( b_prov_s[prov(n)] + Σ_j b_prov_slope_s[prov(n),j] .* X(:, rs_prov_idx(j)) )
%   V_Var(s)  = Var_n( b_var_s[var(n)]  + Σ_j b_var_slope_s[var(n),j]  .* X(:, rs_var_idx(j))  )
%   V_Eps(s)  = mean_n( (nu_s/(nu_s-2)) * sigma_s^2 * exp( z_n(s) ) ), z_n=Z_hs*gamma_s（无异方差则0）
%   share_k(s) = V_k(s) / Σ_k V_k(s) ；对 s 取均值并给 95%CI。

% ---------- 参数 ----------
p = inputParser;
addParameter(p,'dpi',300,@(x)isnumeric(x)&&x>=72);
addParameter(p,'chains',4,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'samples_per_chain',[],@(x)isnumeric(x)&&isscalar(x)&&x>=100);
addParameter(p,'by_column',true,@islogical);     % true: 每目标内归一化；false: 全局
addParameter(p,'show_ci',true,@islogical);
addParameter(p,'palette','mainA30',@ischar);
addParameter(p,'font',10,@(x)isnumeric(x)&&isscalar(x));
parse(p,varargin{:}); A = p.Results;

if nargin<1||isempty(results_root), results_root='results_hbm'; end
cachep = fullfile(results_root,'cache','hbm_cache.mat');
assert(isfile(cachep),'找不到缓存：%s（先运行 prep_hbm_inputs / fit / export）', cachep);
S = load(cachep); cache = S.cache;

% --- 基本量 ---
N  = cache.N;  X = cache.X_fe;  K = size(X,2);
P  = cache.P;  V = cache.V;
prov_id = cache.prov_id(:);
var_id  = cache.var_id(:);
targets = string(fieldnames(cache.Y));
Tg = numel(targets);

% --- 颜色 ---
try, C = smf_colorlib(A.palette);
catch, C = [0.25 0.49 0.79; 0.36 0.69 0.82; 0.85 0.39 0.36; 0.80 0.80 0.80];
end
colFE = pickC(C,1,[0.40 0.40 0.40]);
colProv = pickC(C,2,[0.36 0.69 0.82]);
colVar = pickC(C,3,[0.85 0.39 0.36]);
colEps = [0.86 0.86 0.86];

% --- 存放份额的后验汇总 ---
meanShare = nan(Tg,4);  % [FE Prov Var Eps]
loShare   = nan(Tg,4);
hiShare   = nan(Tg,4);

for ti = 1:Tg
    yname = targets(ti);
    ydir  = fullfile(results_root,'M0_main', yname);
    fdraw = fullfile(ydir,'draws_raw.mat');
    fmeta = fullfile(ydir,'fe_meta.mat');
    assert(isfile(fdraw)&&isfile(fmeta), '缺少 %s 或 %s', fdraw,fmeta);

    L = load(fdraw); draws = L.draws; nRows = height(draws);
    M = load(fmeta); meta = M.meta;

    % 选择保留样本（按链截取）
    Cchain = A.chains;
    if isempty(A.samples_per_chain)
        Skeep = floor(nRows / Cchain);
    else
        Skeep = min(A.samples_per_chain, floor(nRows / Cchain));
    end
    nUse = Cchain * Skeep; rowUse = (1:nUse)';

    % 提取必要列
    alpha = pick_simple(draws, {'alpha'}); alpha = alpha(rowUse); %#ok<NASGU>
    sigma = pick_simple(draws, {'sigma'}); sigma = sigma(rowUse);
    nu    = pick_simple(draws, {'nu'});    if isempty(nu), nu = 7*ones(nUse,1); else, nu = nu(rowUse); end

    [beta,  Kgot] = collect_indexed(draws,'beta');   beta = beta(rowUse,1:min(Kgot,K));
    [bP,    Pgot] = collect_indexed(draws,'b_prov'); bP   = pad_to(bP(rowUse,:), nUse, P);
    [bV,    Vgot] = collect_indexed(draws,'b_var');  bV   = pad_to(bV(rowUse,:), nUse, V);

    % 随机斜率
    JP = numel(meta.rand_slope_prov_idx); JV = numel(meta.rand_slope_var_idx);
    BP_S = zeros(nUse,P,JP);
    BV_S = zeros(nUse,V,JV);
    if JP>0
        for j=1:JP
            tmp = nan(nUse,P);
            for p_i=1:P
                nm = sprintf('b_prov_slope_%d_%d', j, p_i);
                if any(strcmpi(draws.Properties.VariableNames, nm))
                    tmp(:,p_i) = draws.(nm)(rowUse);
                else
                    tmp(:,p_i) = 0;
                end
            end
            BP_S(:,:,j) = tmp;
        end
    end
    if JV>0
        for j=1:JV
            tmp = nan(nUse,V);
            for v_i=1:V
                nm = sprintf('b_var_slope_%d_%d', j, v_i);
                if any(strcmpi(draws.Properties.VariableNames, nm))
                    tmp(:,v_i) = draws.(nm)(rowUse);
                else
                    tmp(:,v_i) = 0;
                end
            end
            BV_S(:,:,j) = tmp;
        end
    end

    % 异方差
    H = numel(meta.heteroskedastic_idx);
    Gamma = [];
    if H>0
        Gamma = nan(nUse,H);
        for j=1:H
            nm = sprintf('gamma_hs_%d', j);
            if any(strcmpi(draws.Properties.VariableNames, nm))
                Gamma(:,j) = draws.(nm)(rowUse);
            else
                Gamma(:,j) = 0;
            end
        end
    end

    % 可用行子样本
    y_tr = cache.Y.(yname);
    ok   = isfinite(y_tr) & all(isfinite(X),2);
    idx_ok = find(ok);
    Xuse   = X(ok,:);
    prov_u = prov_id(ok);
    var_u  = var_id(ok);

    % 预生成选择矩阵
    Nuse = numel(idx_ok);
    SelP = sparse(1:Nuse, prov_u, 1, Nuse, P);
    SelV = sparse(1:Nuse, var_u,  1, Nuse, V);

    rsP = meta.rand_slope_prov_idx(:)';
    rsV = meta.rand_slope_var_idx(:)';
    Zhs = [];
    if H>0, Zhs = Xuse(:, meta.heteroskedastic_idx); end

    % 抽样逐一计算
    shares = nan(nUse,4);
    for s = 1:nUse
        fe_vec = Xuse * beta(s,:)';
        V_FE   = var(fe_vec, 1, 'omitnan');

        prov_int = SelP * bP(s,:)';
        prov_sl  = zeros(Nuse,1);
        for j=1:JP
            bj = BP_S(s,:,j)';                 % P×1
            prov_sl = prov_sl + (Xuse(:,rsP(j)) .* (SelP * bj));
        end
        V_Prov = var(prov_int + prov_sl, 1, 'omitnan');

        var_int = SelV * bV(s,:)';
        var_sl  = zeros(Nuse,1);
        for j=1:JV
            bj = BV_S(s,:,j)';                 % V×1
            var_sl = var_sl + (Xuse(:,rsV(j)) .* (SelV * bj));
        end
        V_Var = var(var_int + var_sl, 1, 'omitnan');

        if H>0, zeta = Zhs * Gamma(s,:)'; else, zeta = 0; end
        fac = max(nu(s)/(nu(s)-2), 0);
        eps_var_obs = fac * (sigma(s)^2) .* exp(zeta);
        V_Eps = mean(eps_var_obs, 'omitnan');

        Vsum = V_FE + V_Prov + V_Var + V_Eps;
        if Vsum<=0 || ~isfinite(Vsum), Vsum = eps; end
        shares(s,:) = [V_FE, V_Prov, V_Var, V_Eps] / Vsum;
    end

    % 后验汇总
    meanShare(ti,:) = mean(shares,1,'omitnan');
    prc = prctile(shares,[2.5 97.5],1);
    loShare(ti,:)   = prc(1,:);
    hiShare(ti,:)   = prc(2,:);
end

% 若按“全局”份额
if ~A.by_column
    meanShare = meanShare ./ max(sum(meanShare,2), eps);
    loShare   = loShare   ./ max(sum(loShare,2),   eps);
    hiShare   = hiShare   ./ max(sum(hiShare,2),   eps);
end

% ---------------- 新增：导出 CSV + 打印 ----------------
outdir = fullfile(results_root,'M0_main','_figs');
if ~isfolder(outdir), mkdir(outdir); end

% 组织列（不使用点式动态赋值，仅矩阵拼接）
% 列顺序：Target, FE_mean/lo/hi, ProvRE_mean/lo/hi, VarRE_mean/lo/hi, Residual_mean/lo/hi
col_block = [ ...
    meanShare(:,1) loShare(:,1) hiShare(:,1) ...
    meanShare(:,2) loShare(:,2) hiShare(:,2) ...
    meanShare(:,3) loShare(:,3) hiShare(:,3) ...
    meanShare(:,4) loShare(:,4) hiShare(:,4) ];

varNames = { ...
    'Target','FE_mean','FE_lo','FE_hi', ...
    'ProvRE_mean','ProvRE_lo','ProvRE_hi', ...
    'VarRE_mean','VarRE_lo','VarRE_hi', ...
    'Residual_mean','Residual_lo','Residual_hi'};

Tcsv = table( cellstr(targets(:)), ...
              col_block(:,1), col_block(:,2), col_block(:,3), ...
              col_block(:,4), col_block(:,5), col_block(:,6), ...
              col_block(:,7), col_block(:,8), col_block(:,9), ...
              col_block(:,10),col_block(:,11),col_block(:,12), ...
              'VariableNames', varNames );

writetable(Tcsv, fullfile(outdir,'Fig1d_variance_shares.csv'));
fprintf('[Fig1d] Variance shares (posterior mean [95%% CrI])\n');
for i=1:Tg
    fprintf('  %-12s | FE %.2f [%.2f,%.2f] | ProvRE %.2f [%.2f,%.2f] | VarRE %.2f [%.2f,%.2f] | Resid %.2f [%.2f,%.2f]\n', ...
        char(targets(i)), ...
        meanShare(i,1), loShare(i,1), hiShare(i,1), ...
        meanShare(i,2), loShare(i,2), hiShare(i,2), ...
        meanShare(i,3), loShare(i,3), hiShare(i,3), ...
        meanShare(i,4), loShare(i,4), hiShare(i,4));
end

% ---------------- 画图（可选中编辑） ----------------
fig = figure('Name','Variance attribution','Color','w','Position',[100 100 1100 720]);
ax  = axes('Parent',fig,'Units','normalized','Position',[0.14 0.15 0.78 0.72]); hold(ax,'on'); safe_pick(ax);

Y = 1:Tg; order = Tg:-1:1;
B = barh(ax, Y, meanShare(order,:), 'stacked');
set(B(1),'FaceColor',colFE,  'EdgeColor','none');
set(B(2),'FaceColor',colProv,'EdgeColor','none');
set(B(3),'FaceColor',colVar, 'EdgeColor','none');
set(B(4),'FaceColor',colEps, 'EdgeColor','none');

set(ax,'YTick',Y,'YTickLabel',targets(order),'TickLabelInterpreter','tex', ...
       'XLim',[0 1], 'Box','off', 'FontSize',A.font*2);
xlabel(ax,'Share of variance (transformed scale)','Interpreter','none'); safe_pick(get(ax,'XLabel'));
ylabel(ax,'Target','Interpreter','none');                                  safe_pick(get(ax,'YLabel'));

if A.show_ci
    for i=1:Tg
        x0 = 0;
        for k=1:4
            m = meanShare(order(i),k);
            lo= loShare(order(i),k);
            hi= hiShare(order(i),k);
            xL = x0 + lo; xH = x0 + hi;
            line(ax,[xL xH],[Y(i) Y(i)],'Color',[0 0 0]*0.45,'LineWidth',0.8,'HitTest','on');
            x0 = x0 + m;
        end
    end
end

L = legend(ax,{'Fixed effects','Province RE','Variant RE','Residual'}, ...
    'Location','southoutside','Orientation','horizontal','Box','off'); safe_pick(L);
ht = title(ax,'Variance attribution by target','FontSize',14,'FontWeight','bold','Interpreter','none'); safe_pick(ht);

if A.dpi>0
    print(fig, fullfile(outdir,'Fig1d_variance_attrib.png'), '-dpng', sprintf('-r%d',A.dpi));
end
end

% ---------------- 工具函数 ----------------
function col = pickC(C,i,fallback)
if size(C,1)>=i, col=C(i,:); else, col=fallback; end
end

function Mout = pad_to(Min,nrow,ncol)
if isempty(Min), Mout=zeros(nrow,ncol); return; end
[r,c]=size(Min); Mout=zeros(nrow,ncol);
Mout(:,1:min(c,ncol)) = Min(:,1:min(c,ncol));
end

function x = pick_simple(T,cand)
x=[]; vn=string(T.Properties.VariableNames);
for i=1:numel(cand)
    hit=find(strcmpi(vn,cand{i}),1);
    if ~isempty(hit), x=T.(vn(hit)); return; end
end
end

function [M,count] = collect_indexed(T,base)
vars=string(T.Properties.VariableNames); ord=nan(size(vars)); keep=false(size(vars));
pat="(?i)"+base+".*?(\d+)(?!.*\d)";
for k=1:numel(vars)
    tok=regexp(vars(k),pat,'tokens','once');
    if ~isempty(tok), ord(k)=str2double(tok{1}); keep(k)=true; end
end
names_all=vars(keep); ord_all=ord(keep);
[ord_sorted,order]=sort(ord_all);
names=names_all(order);
[~,ia]=unique(ord_sorted,'stable'); names=names(ia); ord_sorted=ord_sorted(ia);
count=numel(names);
if count==0, M=zeros(height(T),0); return; end
Sdraw=height(T); M=nan(Sdraw,count);
for j=1:count, M(:,j)=T.(names(j)); end
end

function safe_pick(h)
try
    if isgraphics(h)
        if isprop(h,'HitTest'),       set(h,'HitTest','on'); end
        if isprop(h,'PickableParts'), set(h,'PickableParts','all'); end
    end
catch, end
end
