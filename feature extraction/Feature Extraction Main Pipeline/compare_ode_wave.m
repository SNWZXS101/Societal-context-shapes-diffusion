function compare_ode_wave(matPath, sel, newRootDir)
% compare_ode_wave('waves_dataset.mat', 'Hebei-01', 'D:\...\210waves')
%
% Inputs:
%   matPath    - path to waves_dataset.mat
%   sel        - which wave to compare: 'Hebei-01' or global index (e.g. 123)
%   newRootDir - (optional) new base dir where 210waves are stored.
%                e.g. 'D:\GMDX\...\210waves'

% --- load dataset ---
S = load(matPath, 'Waves'); Waves = S.Waves;

% --- resolve selection ---
if isnumeric(sel) && isscalar(sel)
    k = sel;
elseif ischar(sel) || isstring(sel)
    wid = string(sel);
    % 用 sprintf 拼接 wave_id 再比对
    all_ids = arrayfun(@(W) sprintf('%s-%02d', W.province, W.wave_idx), Waves, 'UniformOutput', false);
    k = find(strcmpi(all_ids, wid), 1);
else
    error('Selector must be numeric index or "Province-XX".');
end
if isempty(k), error('Target wave not found.'); end

W = Waves(k);
wave_id = sprintf('%s-%02d', W.province, W.wave_idx);
fprintf('Comparing %s (variant=%s, beta=%.4f)\n', wave_id, W.variant, W.beta);

% --- determine Excel path ---
if nargin>=3 && ~isempty(newRootDir)
    % 构造新路径：newRootDir\Province\ResultsXX\results\time_series.xlsx
    resultsDir = sprintf('Results%02d', W.wave_idx);
    ts_path = fullfile(newRootDir, W.province, resultsDir, 'results', 'time_series.xlsx');
else
    % 用原始 source_path
    ts_path = W.source_path;
end
if ~isfile(ts_path)
    error('time_series.xlsx not found: %s', ts_path);
end

% --- read original I_ode from Excel ---
Tstab = readtable(ts_path);
vn = lower(string(matlab.lang.makeValidName(Tstab.Properties.VariableNames)));
i_Iode = find(strcmp(vn,'i_ode'),1);
if isempty(i_Iode)
    error('Column I_ode not found in %s', ts_path);
end
I_ode_orig = double(Tstab{:, i_Iode});
I_ode_orig = I_ode_orig(:);

% --- recompute I_ode using integrator ---
gamma = 1/7;
traj = integrate_seir_odeint_style(W.f, W.init, W.beta, gamma);
I_ode_new = traj.I_ode(:);

% --- align length ---
Tmin = min(numel(I_ode_orig), numel(I_ode_new));
I_ode_orig = I_ode_orig(1:Tmin);
I_ode_new  = I_ode_new(1:Tmin);
t0 = (0:Tmin-1)';

% --- plot ---
figure('Color','w');
plot(t0, I_ode_orig, 'o', 'DisplayName','I\_ode (Excel 原始)'); hold on;
plot(t0, I_ode_new, '-', 'LineWidth',1.5, 'DisplayName','I\_ode (回代)');
grid on; xlabel('Day (0-based)'); ylabel('I\_ode');
title(sprintf('%s | %s', wave_id, W.variant));
legend('Location','best');
end
