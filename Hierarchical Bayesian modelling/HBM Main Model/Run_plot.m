% plot_hbm_results('results_root','results_hbm');  % 或者省略参数用默认


% plot_hbm_results('results_root','results_hbm', ...
%     'label_map_json','results_hbm/labels_custom.json', ...
%     'select_covariates', {'gdp_pc','urban population density'}, ...
%     'top_k_coeff', 8);


plot_hbm_results('results_root','results_hbm', 'targets', {'delta_early','kappa_cross','rho_eps','tf1_mean_early','tf2_peak','t_tf2_peak'});
