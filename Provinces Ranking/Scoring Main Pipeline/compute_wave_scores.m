function WaveScores = compute_wave_scores(features_mat, cf_results_path, varargin)
% COMPUTE_WAVE_SCORES  合并 features_aug 与反事实 %Avoided（指定场景），计算逐波 RRI/MEI
% 用法：
%   WaveScores = compute_wave_scores('features_aug.mat', 'cf_results', ...
%       'scenario','TGC_R2ZERO','outdir','scores');
%
% 参数：
%   'scenario'  : 评分采用的反事实场景（默认 'TGC_R2ZERO'）
%   'outdir'    : 输出目录（默认 'scores'）
%   'timestamp' : 是否加时间戳（默认 true）

p = inputParser; p.FunctionName = mfilename;
addParameter(p,'scenario','TGC_STRONG',@(x)ischar(x)||isstring(x));
addParameter(p,'outdir','scores',@(x)ischar(x)||isstring(x));
addParameter(p,'timestamp',true,@(x)islogical(x)||ismember(x,[0 1]));

parse(p,varargin{:});
A = p.Results;

T = load_features_aug(features_mat);
Ta = load_cf_avoided(cf_results_path, A.scenario);

% 合并
K = innerjoin(T, Ta, 'Keys', {'province','wave_id'});

% 需要的列
need = {'delta_early','kappa_cross','tf1_mean_early','t_cross', ...
        'rho_eps','tf2_peak','t_tf2_peak','pct_avoided_main', ...
        'w_tail','flag_multipart'};
check_cols(K, need);

% --- component scores (0–100), 越大越好 ---
n = height(K);
% 趋向越小越好 → 取负号后打分
K.score_delta_early = qscore(-K.delta_early, n);
K.score_kappa       = qscore(-K.kappa_cross, n);
K.score_tf1_early   = qscore(-K.tf1_mean_early, n);
K.score_rho_eps     = qscore( K.rho_eps, n);
K.score_tf2_peak    = qscore(-K.tf2_peak, n);
K.score_t_tf2_peak  = qscore(-K.t_tf2_peak, n);
K.score_pct_avoided = qscore( K.pct_avoided_main, n);

% 惩罚项
lambda_multi = 3 * double(logical(K.flag_multipart));

% “未入负相”的惩罚（t_cross 为 NaN 或非正）
lambda_neg = zeros(n,1);
if ismember('t_cross', K.Properties.VariableNames)
    lambda_neg = 5 * double( isnan(K.t_cross) | (K.t_cross<=0) );
end

% 指数
K.RRI_wave = 0.40*K.score_delta_early + 0.30*K.score_kappa + 0.30*K.score_tf1_early ...
             - lambda_neg - lambda_multi;

K.MEI_wave = 0.30*K.score_rho_eps + 0.10*K.score_tf2_peak + 0.10*K.score_t_tf2_peak ...
             + 0.50*K.score_pct_avoided - lambda_multi;

% 组内排名（分数降序 → 名次升序）
K.rank_RRI_wave = tiedrank(-K.RRI_wave);
K.rank_MEI_wave = tiedrank(-K.MEI_wave);

% 输出精简列（且保留追溯字段）
keep = ["province","wave_id","variant","w_tail","flag_multipart","t_cross", ...
        "delta_early","kappa_cross","tf1_mean_early","rho_eps","tf2_peak","t_tf2_peak","pct_avoided_main", ...
        "score_delta_early","score_kappa","score_tf1_early","score_rho_eps","score_tf2_peak","score_t_tf2_peak","score_pct_avoided", ...
        "RRI_wave","MEI_wave","rank_RRI_wave","rank_MEI_wave"];
keep = keep(ismember(keep, K.Properties.VariableNames));
WaveScores = K(:, keep);

% 落盘
outdir = char(A.outdir); if ~isfolder(outdir), mkdir(outdir); end
if A.timestamp
    mfile = fullfile(outdir, sprintf('wave_scores.mat'));
    cfile = fullfile(outdir, sprintf('wave_scores.csv'));
else
    mfile = fullfile(outdir, 'wave_scores.mat');
    cfile = fullfile(outdir, 'wave_scores.csv');
end
save(mfile, 'WaveScores');
writetable(WaveScores, cfile);
fprintf('[compute_wave_scores] Saved: %s\n', cfile);
end

% ---------- helpers ----------
function check_cols(T, cols)
for i = 1:numel(cols)
    if ~ismember(cols{i}, T.Properties.VariableNames)
        error('缺少列：%s', cols{i});
    end
end
end

function s = qscore(x, n)
% 分位打分：0–100，ties 用平均秩；n=全样本数（210）
if nargin<2, n = numel(x); end
rk = tiedrank(x);  % 最小=1
if n<=1
    s = 50*ones(size(x));
else
    s = 100 * (rk - 1) / (n - 1);
end
end
