



counterfactual_main('features.mat','results_cf', ...
    'auto_T', true, 'T_h', 7, ...
    'neg_floor_ratio',1.0, 'update_mask_mode','adaptive_floor');


S  = load(fullfile('results_cf','timeseries','wave_Anhui-01.mat'));  % 换成你的某一波
TS = S.TS;

fprintf('TGC_STR: Th_used=%g, t0=%g, min(f_TGC)=%.3f, min(f_base)=%.3f\n', ...
    TS.meta.Th_TGC_STRONG, TS.meta.t0_TGC_STRONG, min(TS.f_TGC_STRONG), min(TS.f));



plot_cf_combo('results_cf', 'Anhui-01', ...
    'Scenarios', {'TGC_STRONG','NEG_STRONG'}, ...3
    'SavePNG', false);



% TGC_R2ZERO,TGC_STRONG  TGC_LONG  ADV3


% S = load('features.mat','f_cell','init_cell','Features');
% 
% k    = 1;                                 % 任意一波
% f0   = double(S.f_cell{k}(:));            % 列向量
% init = S.init_cell{k};                    % struct，含 Se0/S0/I0/R0
% beta = S.Features.beta(k);
% if ismember('gamma', S.Features.Properties.VariableNames)
%     gamma = S.Features.gamma(k);
% else
%     gamma = 1/7;
% end
% 
% traj = integrate_seir_odeint_style(f0, init, beta, gamma);  % 这行若能跑通，主程序就能通
% plot(traj.t, traj.I); grid on; title('I(t) baseline');
