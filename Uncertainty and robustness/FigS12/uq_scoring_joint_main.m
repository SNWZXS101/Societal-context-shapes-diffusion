function uq_scoring_and_robustness_main()
% UQ_SCORING_AND_ROBUSTNESS_MAIN
% Scoring/LPI 与 Joint（到 SRS/LPI 层）的不确定性与稳健性：
% - 惩罚稳健性： (lambda_neg, lambda_multi) ∈ {(5,3),(3,2),(7,5)}
% - tail 稳健性： k ∈ {7,10,14} 天窗口（若缺列则现场计算并写回 Features）
% - 权重不确定性：SRS 的 (w_RRI,w_MEI) ~ Dirichlet(0.4,0.6)×scale
%                 LPI 的 (w_delta,w_slope) ~ Dirichlet(0.5,0.5)×scale
% 最终将“惩罚×tail×权重抖动”合并为一个混合抽样 → 有信息量的 rank 分布

%% ---------------- 配置 ----------------
features_mat   = 'features_aug.mat';       % 你真实的 features 文件（会被写回三种 w_tail_k）
cf_results_dir = 'results_cf';
base_outdir    = 'scores';

root_uq  = fullfile('results_hbm_uq');
out_csv  = fullfile(root_uq,'uq_scores');
out_fig  = fullfile(root_uq,'figs_scores_uq');
if ~isfolder(root_uq), mkdir(root_uq); end
if ~isfolder(out_csv), mkdir(out_csv); end
if ~isfolder(out_fig), mkdir(out_fig); end

% 惩罚列表
penalty_list = [5 3; 3 2; 7 5];          % [lambda_neg, lambda_multi]

% Dirichlet 集中程度（越大越贴近名义权重）
alpha_SRS_scale = 150;   % SRS 权重抖动
alpha_LPI_scale = 150;   % LPI 权重抖动
N_weight_draws  = 3000;   % 每个规格下的权重抽样次数
rng(20251029);

% 图形参数
set(0,'defaultTextInterpreter','none');
set(0,'defaultAxesFontSize',11);

%% -------- A. 确保 features 中有 w_tail_7/10/14（没有就现场计算） ----------
add_tail_weights_to_features_if_missing(features_mat);

%% -------- B. 主分析一次，建立省份顺序与基础表 -------------------------
WaveScores0 = compute_wave_scores(features_mat, cf_results_dir, ...
    'scenario','TGC_STRONG', 'outdir',base_outdir, 'timestamp',false);
ProvScores0 = residualize_value_added(WaveScores0, ...
    'use_variant',true, 'weights_col','w_tail', ...
    'w_RRI',0.4,'w_MEI',0.6, ...
    'outdir',base_outdir, 'timestamp',false, 'verbose',true);
prov_list = string(ProvScores0.province);

%% -------- C. 注入 tail 变体（把 7/10/14 并入 WaveScores） ----------------
[WaveScores0, tail_cols] = inject_tail_variants_from_features(WaveScores0, features_mat);
if isempty(tail_cols)
    % 理论上不会，因为上一步已写回；兜底：使用现有 w_tail（视作14天）
    tail_cols = {'w_tail'};
end
% 标准化 tail 标签（用于图注/文件名）
tail_label = @(colname) strrep(strrep(char(colname),'w_tail','tail='),'_','');

%% -------- D. 构造“规格-权重混合抽样” -----------------------------------
% 对每个 (惩罚, tail) 规格：
%   1) 复算 RRI_wave/MEI_wave；
%   2) 以该 tail 权重聚合，得到省级 VA（RRI_VA/MEI_VA）；
%   3) 对 SRS 权重做 Dirichlet 抖动得到 SRS 抽样；
%   4) 基于 wave 轨迹 + LPI 组件，对 LPI 权重抖动得到 LPI 抽样；
% 然后把所有规格的抽样在“行维”串接，形成一个大抽样矩阵。

SRS_draws_all = [];
LPI_draws_all = [];
spec_names     = {};   % 仅用于跨规格热图

ip = 1;
while ip<=size(penalty_list,1)
    lam_neg   = penalty_list(ip,1);
    lam_multi = penalty_list(ip,2);

    it = 1;
    while it<=numel(tail_cols)
        tcol = tail_cols{it};
        tag  = sprintf('pen=(%d,%d) %s',lam_neg,lam_multi,tail_label(tcol));

        % (1) 惩罚复算
        WS_pen = recompute_wave_scores_with_penalty(WaveScores0, lam_neg, lam_multi);

        % (2) tail 权重覆盖
        if ~ismember(tcol, WS_pen.Properties.VariableNames)
            warning('未在 WaveScores 中找到列 %s，跳过该 tail 规格。', tcol);
            it = it + 1; continue;
        end
        WS_pen.w_tail = WS_pen.(tcol);

        % (3) SRS 抖动抽样
        [~, SRS_sum_tbl, ~, ~, SRS_draws_spec] = srs_with_dirichlet_draws_spec( ...
            WS_pen, prov_list, alpha_SRS_scale, N_weight_draws);
        % (4) LPI 抖动抽样
        [~, LPI_sum_tbl, ~, ~, LPI_draws_spec] = lpi_with_dirichlet_draws_spec( ...
            WS_pen, prov_list, features_mat, alpha_LPI_scale, N_weight_draws);

        SRS_draws_all = [SRS_draws_all; SRS_draws_spec]; %#ok<AGROW>
        LPI_draws_all = [LPI_draws_all; LPI_draws_spec]; %#ok<AGROW>
        spec_names{end+1} = tag; %#ok<AGROW>

        % 落盘每个规格的 summary（可选）
        writetable(SRS_sum_tbl, fullfile(out_csv, ['SRS_' sanitize(tag) '_summary.csv']));
        writetable(LPI_sum_tbl, fullfile(out_csv, ['LPI_' sanitize(tag) '_summary.csv']));

        it = it + 1;
    end
    ip = ip + 1;
end

% 抽样矩阵 → 汇总与排名（混合后的“总体不确定性”）
[~, SRS_sum, SRS_rank, SRS_ridge] = summarize_draws('SRS', prov_list, SRS_draws_all);
[~, LPI_sum, LPI_rank, LPI_ridge] = summarize_draws('LPI', prov_list, LPI_draws_all);
writetable(SRS_sum,  fullfile(out_csv,'SRS_mixture_summary.csv'));
writetable(SRS_rank, fullfile(out_csv,'SRS_mixture_ranks.csv'));
writetable(SRS_ridge,fullfile(out_csv,'SRS_mixture_rankridge.csv'));
writetable(LPI_sum,  fullfile(out_csv,'LPI_mixture_summary.csv'));
writetable(LPI_rank, fullfile(out_csv,'LPI_mixture_ranks.csv'));
writetable(LPI_ridge,fullfile(out_csv,'LPI_mixture_rankridge.csv'));

%% -------- E. 可视化（新的三类图 + 跨规格均值热图） ----------------------
% 1) 横向误差条（Top/Bottom 10）
plot_barh_ci(SRS_draws_all, prov_list, 'SRS | mixture (penalties × tails × weights)', out_fig);
plot_barh_ci(LPI_draws_all, prov_list, 'LPI | mixture (penalties × tails × weights)', out_fig);

% 2) Rank 概率热图（非退化）
plot_rank_heatmap(SRS_draws_all, prov_list, 'SRS | mixture rank probability', out_fig);
plot_rank_heatmap(LPI_draws_all, prov_list, 'LPI | mixture rank probability', out_fig);

% 3) Rank ridge（山峦）图
plot_rank_ridge(SRS_draws_all, prov_list, 'SRS | mixture rank ridge', out_fig);
plot_rank_ridge(LPI_draws_all, prov_list, 'LPI | mixture rank ridge', out_fig);

% 4) 跨规格均值热图（行=省，列=规格；值=该规格下的均值）
plot_crossspec_mean_heatmap(SRS_draws_all, prov_list, spec_names, 'SRS | cross-spec means', out_fig);
plot_crossspec_mean_heatmap(LPI_draws_all, prov_list, spec_names, 'LPI | cross-spec means', out_fig);

disp('== UQ 完成，见 results_hbm_uq/uq_scores 与 results_hbm_uq/figs_scores_uq ==');
end

%% ===================== 工具函数 ===========================

function add_tail_weights_to_features_if_missing(features_mat)
% 若 Features 中缺 w_tail_7/10/14 → 依据 T_tail 公式现场计算并写回
if ~isfile(features_mat), error('找不到 %s', features_mat); end
S = load(features_mat);
Tfeat = [];
names = {'Features','Tfeat','features','T'};
tname = '';
i=1;
while i<=numel(names)
    if isfield(S,names{i}), Tfeat = S.(names{i}); tname = names{i}; break; end
    i=i+1;
end
if ~istable(Tfeat), error('features_mat 中未找到表（Features/Tfeat/…）'); end

cols = string(Tfeat.Properties.VariableNames);
has7  = any(strcmpi(cols,'w_tail_7') | strcmpi(cols,'w_tail7'));
has10 = any(strcmpi(cols,'w_tail_10')| strcmpi(cols,'w_tail10'));
has14 = any(strcmpi(cols,'w_tail_14')| strcmpi(cols,'w_tail14'));

if has7 && has10 && has14
    fprintf('[features] 已存在 w_tail_7/10/14，跳过计算。\n'); return;
end

% 宽松找列名
cT   = pick_col(cols, {'T','t_len','tL','t_l','t_total'});
ctf  = pick_col(cols, {'t_fmax','tf2_peak','t_tf2_peak','t_f_max'});
cinc = pick_col(cols, {'t_inc_slope_max','t_inc_slope','t_incmax','t_inc_max'});
if any(strcmp('',{cT,ctf,cinc}))
    error('Features 缺 T / t_fmax / t_inc_slope_max（或同义列）。');
end

T  = double(Tfeat.(cT));
tf = double(Tfeat.(ctf));
ti = double(Tfeat.(cinc));
T_tail = T - max(tf,ti);
T_tail(~isfinite(T_tail)) = NaN;
T_tail(T_tail<0) = 0;

if ~has7,  Tfeat.w_tail_7  = min(1, T_tail/7);  end
if ~has10, Tfeat.w_tail_10 = min(1, T_tail/10); end
if ~has14, Tfeat.w_tail_14 = min(1, T_tail/14); end

S.(tname) = Tfeat;
save(features_mat,'-struct','S');
fprintf('[features] 已写入 w_tail_7/10/14 到 %s（表：%s）。\n', features_mat, tname);
end

function [WS2, tail_cols] = inject_tail_variants_from_features(WaveScores, features_mat)
% 将不同写法的 tail 列并入 WaveScores，标准化列名为 w_tail7/10/14
WS2 = WaveScores; tail_cols = {};
if ~isfile(features_mat), return; end
S = load(features_mat);
Tfeat = [];
names={'Features','Tfeat','features','T'}; tname='';
i=1;while i<=numel(names)
    if isfield(S,names{i}), Tfeat=S.(names{i}); tname=names{i}; break; end
    i=i+1;
end
if ~istable(Tfeat), return; end

K = {'province','wave_id'};
for i=1:numel(K)
    if ismember(K{i},Tfeat.Properties.VariableNames)
        if ~isstring(Tfeat.(K{i})), Tfeat.(K{i}) = string(Tfeat.(K{i})); end
    end
    if ismember(K{i},WS2.Properties.VariableNames)
        if ~isstring(WS2.(K{i})), WS2.(K{i}) = string(WS2.(K{i})); end
    end
end

cols = string(Tfeat.Properties.VariableNames);
tok = regexp(lower(cols), '^w[_]?tail[_]?(\d+)$','tokens','once');
want = {'7','10','14'}; std = {'w_tail7','w_tail10','w_tail14'};
src = strings(1,3);
j=1;while j<=numel(cols)
    t = tok{j};
    if ~isempty(t)
        d = string(t{1});
        k = find(strcmp(d,want),1);
        if ~isempty(k) && src(k)==""
            src(k) = cols(j);
        end
    end
    j=j+1;
end

k=1;while k<=3
    if src(k)~=""
        Tuse = Tfeat(:,[K, cellstr(src(k))]);
        Tuse.Properties.VariableNames(end) = {std{k}};
        WS2 = outerjoin(WS2, Tuse, 'Keys',K,'MergeKeys',true);
        tail_cols{end+1} = std{k}; %#ok<AGROW>
    end
    k=k+1;
end

if ~isempty(tail_cols)
    fprintf('[inject] 检测到 tail 变体：%s\n', strjoin(tail_cols, ', '));
else
    fprintf('[inject] 未识别到 w_tail7/10/14（Features:%s）。\n', tname);
end
end

function col = pick_col(cols, aliases)
col = '';
i=1;while i<=numel(aliases)
    m = strcmpi(cols, aliases{i});
    if any(m), col = cols(find(m,1)); return; end
    i=i+1;
end
end

function S = sanitize(tag)
S = regexprep(tag,'[^\w\- =(),]','_');
S = strrep(S,'—','-');
S = strrep(S,' ','_');
end

function WS2 = recompute_wave_scores_with_penalty(WaveScores, lambda_neg, lambda_multi)
% 按给定惩罚复算 RRI_wave/MEI_wave
K = WaveScores;
need = {'score_delta_early','score_kappa','score_tf1_early', ...
        'score_rho_eps','score_tf2_peak','score_t_tf2_peak','score_pct_avoided', ...
        't_cross','flag_multipart'};
for i=1:numel(need)
    if ~ismember(need{i}, K.Properties.VariableNames)
        error('WaveScores 缺列：%s', need{i});
    end
end
n = height(K);
lam_multi = lambda_multi * double(logical(K.flag_multipart));
lam_neg = lambda_neg * double( isnan(K.t_cross) | (K.t_cross<=0) );

K.RRI_wave = 0.40*K.score_delta_early + 0.30*K.score_kappa + 0.30*K.score_tf1_early ...
             - lam_neg - lam_multi;
K.MEI_wave = 0.30*K.score_rho_eps + 0.10*K.score_tf2_peak + 0.10*K.score_t_tf2_peak ...
             + 0.50*K.score_pct_avoided - lam_multi;
K.rank_RRI_wave = tiedrank(-K.RRI_wave);
K.rank_MEI_wave = tiedrank(-K.MEI_wave);
WS2 = K;
end

function [draws_tbl, sum_tbl, rank_tbl, ridge_tbl, Dmat] = srs_with_dirichlet_draws_spec(WaveScores, prov_list, alpha_scale, N)
% 给定一个 (惩罚,tail) 规格，生成 SRS 的 N 次权重抖动抽样
PS = residualize_value_added(WaveScores, 'use_variant',true, ...
    'weights_col','w_tail','w_RRI',0.4,'w_MEI',0.6, ...
    'outdir','scores','timestamp',false,'verbose',false);
PS = sortrows(PS,'province');
[~, ia] = ismember(string(prov_list), string(PS.province));
P = numel(prov_list);
R = nan(P,1); E = nan(P,1); mk = ia>0;
R(mk) = PS.RRI_VA(ia(mk));
E(mk) = PS.MEI_VA(ia(mk));

alpha = alpha_scale * [0.4 0.6];
W = dirichlet_samples(alpha, N);   % N×2
D = N;
SRS = nan(D,P);
i=1;while i<=D
    wR=W(i,1); wE=W(i,2);
    s = (wR*R + wE*E) / max(1e-12,(wR+wE));
    SRS(i,:) = s';
    i=i+1;
end
[draws_tbl, sum_tbl, rank_tbl, ridge_tbl] = summarize_draws('SRS', prov_list, SRS);
Dmat = SRS;
end

function [draws_tbl, sum_tbl, rank_tbl, ridge_tbl, Dmat] = lpi_with_dirichlet_draws_spec(WaveScores, prov_list, features_mat, alpha_scale, N)
% 基于 WaveScores（含该规格的 w_tail），先得到 LPI 组件，再做 N 次权重抖动
PS = residualize_value_added(WaveScores, 'use_variant',true, ...
    'weights_col','w_tail','w_RRI',0.4,'w_MEI',0.6, ...
    'outdir','scores','timestamp',false,'verbose',false);
LP = compute_lpi(WaveScores, PS, ...
    'features_mat',features_mat,'delta_ceiling',true, ...
    'slope_method','theilsen','w_delta',0.5,'w_slope',0.5, ...
    'outdir','scores','timestamp',false,'verbose',false);
LP = sortrows(LP,'province');
[~, ia] = ismember(string(prov_list), string(LP.province));
P = numel(prov_list);
cdel = nan(P,1); cslp = nan(P,1); mk = ia>0;
cdel(mk) = LP.comp_delta(ia(mk));
cslp(mk) = LP.comp_slope(ia(mk));

alpha = alpha_scale * [0.5 0.5];
W = dirichlet_samples(alpha, N);
D = N;
LPI = nan(D,P);
i=1;while i<=D
    w1=W(i,1); w2=W(i,2);
    y = w1*cdel + w2*cslp;
    LPI(i,:) = y';
    i=i+1;
end
[draws_tbl, sum_tbl, rank_tbl, ridge_tbl] = summarize_draws('LPI', prov_list, LPI);
Dmat = LPI;
end

function [draws_tbl, sum_tbl, rank_tbl, ridge_tbl] = summarize_draws(~, prov, Dmat)
D = size(Dmat,1); P = size(Dmat,2);
mu = mean(Dmat,1,'omitnan');
lo = quantile(Dmat,0.025,1);
hi = quantile(Dmat,0.975,1);

% 稳健排名（NaN 设为最差）
R = zeros(D,P);
i=1;while i<=D
    v = Dmat(i,:);
    m = isfinite(v);
    rk = P*ones(1,P);
    if any(m)
        [~,o]=sort(v(m),'descend'); rr=1:sum(m);
        idx=find(m); rk(idx(o)) = rr;
    end
    R(i,:)=rk; i=i+1;
end
rmean = mean(R,1,'omitnan');
rlo   = quantile(R,0.025,1);
rhi   = quantile(R,0.975,1);
topk  = min(10,P);
topk_prob = mean(R<=topk,1,'omitnan');

prov_col = string(prov(:));
draws_tbl = array2table(Dmat, 'VariableNames', cellstr(prov_col));
sum_tbl = table(prov_col, mu(:), lo(:), hi(:), topk_prob(:), ...
    'VariableNames', {'Province','Mean','CrI_lo','CrI_hi','Top10Prob'});
rank_tbl = table(prov_col, rmean(:), rlo(:), rhi(:), ...
    'VariableNames', {'Province','RankMean','RankCrI_lo','RankCrI_hi'});

% rank 概率
RR = zeros(P,P);
j=1;while j<=P
    rr = R(:,j);
    h = histcounts(rr, 0.5:1:(P+0.5));
    RR(j,:) = h / max(1,sum(isfinite(rr)));
    j=j+1;
end
Province = strings(P*P,1); Rank=zeros(P*P,1); Prob=zeros(P*P,1);
cnt=1; a=1;while a<=P
    b=1;while b<=P
        Province(cnt)=prov_col(a); Rank(cnt)=b; Prob(cnt)=RR(a,b);
        cnt=cnt+1; b=b+1;
    end
    a=a+1;
end
ridge_tbl = table(Province, Rank, Prob);
end

function W = dirichlet_samples(alpha, N)
K = numel(alpha);
if exist('gamrnd','file')==2
    G = gamrnd(repmat(alpha,[N,1]),1.0);
elseif exist('randg','file')==2
    G = randg(repmat(alpha,[N,1]));
else
    % 兜底近似
    m = alpha(:)'/sum(alpha);
    G = max(1e-6, repmat(m,[N,1]) + 0.05*randn(N,K));
end
S = sum(G,2);
W = G ./ S;
end

%% -------------------- 绘图 -----------------------------

function plot_barh_ci(Dmat, prov_list, tag, out_dir)
mu = mean(Dmat,1,'omitnan'); [~,ord]=sort(mu,'descend');
X  = Dmat(:,ord); prov2 = prov_list(ord);
lo = quantile(X,0.025,1); hi = quantile(X,0.975,1);
P  = numel(prov2); k = min(10, floor(P/3));
idx = [1:k, (P-k+1):P]; prov3 = prov2(idx);
mu3 = mu(ord(idx)); lo3 = lo(idx); hi3 = hi(idx);

fig = figure('Position',[80 80 950 120+24*numel(idx)],'Renderer','painters');
hold on;
y = 1:numel(idx);
i=1;while i<=numel(idx)
    plot([lo3(i) hi3(i)],[y(i) y(i)],'-','LineWidth',1.6);
    i=i+1;
end
plot(mu3, y, 'o','LineWidth',1.6,'MarkerSize',5);
set(gca,'YTick',y,'YTickLabel',cellstr(prov3));
xlabel('Estimate'); ylabel('Province (top/bottom 10 by mean)');
title([tag ' - horizontal CI bars']);
grid on; box on; drawnow;
export_figs(fig, fullfile(out_dir,[sanitize(tag) '_barh_CI']));
close(fig);
end

function plot_rank_heatmap(Dmat, prov_list, tag, out_dir)
[D,P] = size(Dmat);
% rank 矩阵
R = zeros(D,P);
i=1;while i<=D
    v=Dmat(i,:); m=isfinite(v); rk=P*ones(1,P);
    if any(m)
        [~,o]=sort(v(m),'descend'); rr=1:sum(m); idx=find(m); rk(idx(o))=rr;
    end
    R(i,:)=rk; i=i+1;
end
Prob = zeros(P,P);
j=1;while j<=P
    h = histcounts(R(:,j),0.5:1:(P+0.5));
    Prob(j,:) = h / max(1,sum(isfinite(R(:,j))));
    j=j+1;
end

% 只显示“概率>0.02”的色阶，避免“全蓝”
Prob(Prob<0.02) = 0;

[~,ord]=sort(mean(Dmat,1,'omitnan'),'descend');
prov2 = prov_list(ord);
Prob  = Prob(ord,:);

fig = figure('Position',[80 80 1100 140+20*P],'Renderer','painters');
imagesc(Prob);
set(gca,'YTick',1:P,'YTickLabel',cellstr(prov2),'XTick',1:P,'XTickLabel',1:P);
xlabel('Rank'); ylabel('Province');
colormap(parula); colorbar;
title([tag]);
drawnow;
export_figs(fig, fullfile(out_dir,[sanitize(tag) '_heatmap_rankprob']));
close(fig);
end

function plot_rank_ridge(Dmat, prov_list, tag, out_dir)
[D,P] = size(Dmat);
% rank 概率
R = zeros(D,P);
i=1;while i<=D
    v=Dmat(i,:); m=isfinite(v); rk=P*ones(1,P);
    if any(m)
        [~,o]=sort(v(m),'descend'); rr=1:sum(m); idx=find(m); rk(idx(o))=rr;
    end
    R(i,:)=rk; i=i+1;
end
Prob = zeros(P,P);
j=1;while j<=P
    h = histcounts(R(:,j),0.5:1:(P+0.5));
    Prob(j,:) = h / max(1,sum(isfinite(R(:,j))));
    j=j+1;
end
[~,ord]=sort(mean(Dmat,1,'omitnan'),'descend');
prov2 = prov_list(ord); Prob = Prob(ord,:);

% 画“山峦”：每省一条填充曲线
fig = figure('Position',[80 80 1100 140+18*P],'Renderer','painters');
hold on;
yoff = 0; gap = 1.1; scale = P;  % 纵向间距 & 放缩
p=1;while p<=P
    y = yoff + Prob(p,:)*scale;
    x = 1:P;
    % 填充
    fill([x, fliplr(x)], [yoff+zeros(1,P), fliplr(y)], [0.7 0.7 0.9], ...
        'EdgeColor','k','LineWidth',0.6,'FaceAlpha',0.6);
    text(P+0.5, yoff+0.2*scale, char(prov2(p)), 'HorizontalAlignment','left');
    yoff = yoff + gap*scale;
    p=p+1;
end
xlim([1 P]); ylim([0 yoff]); xlabel('Rank'); yticks([]); box on;
title([tag]);
drawnow;
export_figs(fig, fullfile(out_dir,[sanitize(tag) '_rank_ridge']));
close(fig);
end

function plot_crossspec_mean_heatmap(Dmat_all, prov_list, spec_names, tag, out_dir)
% 将“整体抽样”按规格分块求均值，形成 行=省 / 列=规格 的热图
% 假设每个规格的抽样行数相等（N_weight_draws）
if isempty(spec_names), return; end
S = numel(spec_names);
N_each = size(Dmat_all,1)/S;
if abs(N_each - round(N_each)) > 1e-8
    % 若不整除，退化为整体均值图
    mu = mean(Dmat_all,1,'omitnan')';
    fig = figure('Position',[80 80 900 120+20*numel(prov_list)],'Renderer','painters');
    imagesc(mu'); colormap(parula); colorbar;
    set(gca,'XTick',1,'XTickLabel',{'mean'},'YTick',1:numel(prov_list),'YTickLabel',cellstr(prov_list));
    title([tag ' (overall mean)']); drawnow;
    export_figs(fig, fullfile(out_dir,[sanitize(tag) '_heatmap_overallmean']));
    close(fig);
    return;
end

P = numel(prov_list);
M = zeros(P,S);
s=1;row_start=1;
while s<=S
    rows = row_start:(row_start+N_each-1);
    M(:,s) = mean(Dmat_all(rows,:),1,'omitnan')';
    row_start = row_start + N_each;
    s=s+1;
end

% 以整体均值排序行
mu_all = mean(Dmat_all,1,'omitnan');
[~,ord] = sort(mu_all,'descend');
prov2 = prov_list(ord);
M = M(ord,:);

fig = figure('Position',[80 80 1200 140+20*P],'Renderer','painters');
imagesc(M);
set(gca,'YTick',1:P,'YTickLabel',cellstr(prov2),'XTick',1:S,'XTickLabel',spec_names);
xtickangle(45);
colormap(parula); colorbar;
title([tag]); drawnow;
export_figs(fig, fullfile(out_dir,[sanitize(tag) '_heatmap_crossspec']));
close(fig);
end

function export_figs(fig, base)
% 只保存为可编辑的 .fig
% 使用 'compact' 减小文件体积（R2020a+ 支持；旧版可去掉该参数）
try
    savefig(fig, [base '.fig'], 'compact');
catch
    savefig(fig, [base '.fig']);
end
saveas(fig, [base '.png']);
end
