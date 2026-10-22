function hbm_uq_driver()
clc; fprintf('[HBM-UQ] start\n');

%% 0) 预处理一次（公共 cache）
prep_out = 'results_hbm';
prep_hbm_inputs('out_dir', prep_out, ...
  'keep_targets',{'delta_early','kappa_cross','rho_eps','tf1_mean_early','tf2_peak','t_tf2_peak'}, ...
  'log_targets',  {'delta_early','rho_eps','tf1_mean_early','tf2_peak','t_tf2_peak'}, ...
  'abslog_targets',{'kappa_cross'}, ...
  'hard_drop_nonpos_targets',{'delta_early'}, ...
  'align_xi',false);

% prep_hbm_inputs('out_dir','results_hbm', ...
%   'keep_targets',{'delta_early','kappa_cross','rho_eps','tf1_mean_early','tf2_peak','t_tf2_peak'}, ...
%   'log_targets',  {'delta_early','rho_eps','tf1_mean_early','tf2_peak','t_tf2_peak'}, ...
%   'abslog_targets',{'kappa_cross'}, ...
%   'hard_drop_nonpos_targets',{'delta_early'}, ...  % 仍旧只硬删 delta_early
%   'log_nonpos_to_nan',true, ...
%   'log_warn_on_nonpos',false, ...   % ← 关闭这条提示
%   'align_xi',false);



S = load(fullfile(prep_out,'cache','hbm_cache.mat')); C = S.cache;
targets = fieldnames(C.Y);

%% 1) Define specs（你的版本，已保留）
BASE = struct( ...
  'name','base', ...
  'beta_sd',1, 'alpha_sd',3, 'ranef_scale',3, ...
  'rs_idx',[7 12], 'rs_var_idx',[7 12], 'hs_idx',[8]);
vec0 = @() zeros(0,1);               % 0×1 向量，满足 isnumeric & isvector

specs(1) = BASE;
specs(2) = makeSpec(BASE,'priorsx0p5', 'beta_sd',0.5,'alpha_sd',1.5,'ranef_scale',1.5);
specs(3) = makeSpec(BASE,'priorsx2',   'beta_sd',2.0,'alpha_sd',6.0,'ranef_scale',6.0);
specs(4) = makeSpec(BASE,'no_RS',      'rs_idx',vec0(),'rs_var_idx',vec0());
specs(5) = makeSpec(BASE,'no_HS',      'hs_idx',vec0());
specs(6) = makeSpec(BASE,'no_RS_HS',   'rs_idx',vec0(),'rs_var_idx',vec0(),'hs_idx',vec0());

spec_names = {specs.name};

%% （你漏掉的这两句 ↓↓↓）
root_uq  = 'results_hbm_uq';                      % UQ 根目录
if ~isfolder(root_uq), mkdir(root_uq); end
src_cache = fullfile(prep_out,'cache','hbm_cache.mat');   % 公共 cache 源

%% 2) 逐规格拟合 + 导出
for s = 1:numel(specs)
    cfg = specs(s);
    out_root = fullfile(root_uq, cfg.name);
    if ~isfolder(out_root), mkdir(out_root); end

    % 确保每个规格目录下有 cache/hbm_cache.mat
    cache_dir = fullfile(out_root,'cache');
    if ~isfolder(cache_dir), mkdir(cache_dir); end
    dst_cache = fullfile(cache_dir,'hbm_cache.mat');
    if ~isfile(dst_cache), copyfile(src_cache, dst_cache); end

    fprintf('[HBM-UQ] Fitting spec: %s\n', cfg.name);
    fit_hbm_main_m0_matlab( ...
        'results_root', out_root, ...
        'num_warmup', 4000, 'num_samples', 4000, 'chains', 4, 'seed', 42, ...
        'beta_sd', cfg.beta_sd, 'alpha_sd', cfg.alpha_sd, 'ranef_scale', cfg.ranef_scale, ...
        'learn_nu', false, 'nu_step', 0.10, 'init_nu', 30, 'init_jitter', 0.3, ...
        'rand_slope_prov_idx', cfg.rs_idx, ...
        'rand_slope_var_idx',  cfg.rs_var_idx, ...
        'heteroskedastic_idx', cfg.hs_idx);

    export_hbm_outputs('results_root', out_root, 'n_ppc_draws', 1000, 'chains', 4);
end

%% 3) 汇总并弹窗三张图（诊断 & PPC）
nS = numel(specs); nT = numel(targets);
Rhat_med = nan(nT,nS); ESS_min = nan(nT,nS);
PPC_rho  = nan(nT,nS); PPC_cover = nan(nT,nS);

for s = 1:nS
    out_root = fullfile(root_uq, specs(s).name);
    for ti = 1:nT
        yname = targets{ti}; ydir = fullfile(out_root,'M0_main',yname);

        fsum = fullfile(ydir,'summary.csv');
        if isfile(fsum)
            Sx = readtable(fsum);
            rhat = Sx.rhat; ess = Sx.ess;
            rhat = rhat(isfinite(rhat) & rhat>0);
            ess  = ess(isfinite(ess)  & ess>0);
            if ~isempty(rhat), Rhat_med(ti,s) = median(rhat); end
            if ~isempty(ess),  ESS_min(ti,s)  = min(ess);     end
        end

        fppc = fullfile(ydir,'ppc.csv');
        if isfile(fppc)
            P = readtable(fppc);
            ok  = isfinite(P.y_obs) & isfinite(P.yhat_p50);
            if any(ok),  PPC_rho(ti,s)   = corr(P.y_obs(ok), P.yhat_p50(ok), 'Type','Spearman','Rows','pairwise'); end
            okc = isfinite(P.y_obs) & isfinite(P.yhat_p10) & isfinite(P.yhat_p90);
            if any(okc), PPC_cover(ti,s) = mean(P.y_obs(okc) >= P.yhat_p10(okc) & P.y_obs(okc) <= P.yhat_p90(okc)); end
        end
    end
end

tlabels = strrep(strrep(string(targets),'_','\_'),'delta','\Delta');
tlabels = replace(tlabels,"rho",'\rho'); tlabels = replace(tlabels,"kappa",'\kappa');

% ==== Fig 1: split-Rhat（格内标 min ESS）====
figure('Color','w','Position',[100 100 1100 600]);
imagesc(Rhat_med); axis tight; colormap(flipud(parula)); colorbar;
title('Convergence: split-\hat{R} (median) — numbers show min ESS', ...
      'FontSize',14, 'Interpreter','tex');
set(gca,'XTick',1:nS,'XTickLabel',spec_names,'XTickLabelRotation',35, ...
        'YTick',1:nT,'YTickLabel',tlabels,'FontSize',12,'TickLabelInterpreter','tex');
caxis([0.98 1.02]); hold on;
for i=1:nT, for j=1:nS
    if isfinite(ESS_min(i,j))
        text(j,i,sprintf('%.0f', ESS_min(i,j)), ...
             'HorizontalAlignment','center','VerticalAlignment','middle', ...
             'Color','w','FontWeight','bold','FontSize',11);
    end
end, end
hold off;

% ==== Fig 2: PPC 相关 ====
figure('Color','w','Position',[100 100 1100 600]);
imagesc(PPC_rho); axis tight; colormap(parula); colorbar;
title('Posterior predictive rank correlation  \rho_s(y_{obs}, \hat{y}_{50})', ...
      'FontSize',14,'Interpreter','tex');
set(gca,'XTick',1:nS,'XTickLabel',spec_names,'XTickLabelRotation',35, ...
        'YTick',1:nT,'YTickLabel',tlabels,'FontSize',12,'TickLabelInterpreter','tex');
caxis([0.6 1.0]);
for i=1:nT, for j=1:nS
    if isfinite(PPC_rho(i,j))
        text(j,i,sprintf('%.2f',PPC_rho(i,j)), ...
             'HorizontalAlignment','center','VerticalAlignment','middle', ...
             'Color','k','FontWeight','bold','FontSize',11);
    end
end, end

% ==== Fig 3: PPC 覆盖率（10–90%）====
figure('Color','w','Position',[100 100 1100 600]);
imagesc(PPC_cover); axis tight; colormap(parula); colorbar;
title('Posterior predictive 10–90% coverage', 'FontSize',14, 'Interpreter','tex');
set(gca,'XTick',1:nS,'XTickLabel',spec_names,'XTickLabelRotation',35, ...
        'YTick',1:nT,'YTickLabel',tlabels,'FontSize',12,'TickLabelInterpreter','tex');
caxis([0.6 1.0]);
for i=1:nT, for j=1:nS
    if isfinite(PPC_cover(i,j))
        text(j,i,sprintf('%.2f',PPC_cover(i,j)), ...
             'HorizontalAlignment','center','VerticalAlignment','middle', ...
             'Color','k','FontWeight','bold','FontSize',11);
    end
end, end

end

function S2 = makeSpec(S, name, varargin)
S2 = S; S2.name = char(name);
if ~isempty(varargin)
    assert(mod(numel(varargin),2)==0,'key/value 必须成对');
    for k = 1:2:numel(varargin), S2.(varargin{k}) = varargin{k+1}; end
end
end
