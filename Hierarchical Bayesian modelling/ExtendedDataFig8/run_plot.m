
out_dir = 'results_hbm/figs_burden';
if ~isfolder(out_dir), mkdir(out_dir); end

% Panel A: I_max（蓝）
burden_forest_panel('results_hbm','I_max', 'I_{max} (NB2-GLMM)', ...
    'Kshow',25,'Color',[0.12 0.45 0.75], ...
    'Out', fullfile(out_dir,'PanelA_Imax'));

% Panel B: C_total（绿）
burden_forest_panel('results_hbm','C_total', 'C_{total} (NB2-GLMM)', ...
    'Kshow',25,'Color',[0.10 0.60 0.45], ...
    'Out', fullfile(out_dir,'PanelB_Ctotal'));


out_dir = 'results_hbm/figs_burden';
if ~isfolder(out_dir), mkdir(out_dir); end
% Panel C: 以 I_max 为例；如果想画 C_total，再调用一次改 target_dir 和标题
burden_ppc_panel('results_hbm','I_max', 'Posterior predictive: I_{max}', ...
   'Colors',[0.82 0.33 0.30; 0.95 0.66 0.12; 0.12 0.47 0.71], ...
   'Out', fullfile(out_dir,'PanelC_Imax'));
burden_ppc_panel('results_hbm','C_total', 'Posterior predictive: C_{total}', ...
   'Colors',[0.82 0.33 0.30; 0.95 0.66 0.12; 0.12 0.47 0.71], ...
   'Out', fullfile(out_dir,'PanelC_Ctotal'));


visualize_burden_ED('results_root','results_hbm');



% 
% burden_write_missing_datajson('results_hbm','I_max');
% burden_write_missing_datajson('results_hbm','C_total');
