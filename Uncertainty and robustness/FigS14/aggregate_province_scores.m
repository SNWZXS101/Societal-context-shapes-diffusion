function ProvAgg = aggregate_province_scores(WaveScores, varargin)
% AGGREGATE_PROVINCE_SCORES
% 省内加权聚合 RRI/MEI（含 NaN 策略与日志），输出 ProvAgg 并落盘。
%
% 用法：
%   ProvAgg = aggregate_province_scores(WaveScores, ...
%       'nan_mode','exclude', ...      % 'exclude' | 'prov_mean'
%       'outdir','scores', ...
%       'timestamp',true, ...
%       'verbose',true);

p = inputParser; p.FunctionName = mfilename;
addParameter(p,'nan_mode','exclude',@(s)any(strcmpi(s,{'exclude','prov_mean'})));
addParameter(p,'outdir','scores',@(x)ischar(x)||isstring(x));
addParameter(p,'timestamp',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'verbose',true,@(x)islogical(x)||ismember(x,[0 1]));
parse(p,varargin{:});
A = p.Results;

% ---- 基本检查 ----
need = {'province','w_tail','RRI_wave','MEI_wave'};
for i = 1:numel(need)
    if ~ismember(need{i}, WaveScores.Properties.VariableNames)
        error('WaveScores 缺少列：%s', need{i});
    end
end

% ---- 统一类型 ----
if ~isstring(WaveScores.province), WaveScores.province = string(WaveScores.province); end

% ---- 聚合 ----
G = findgroups(WaveScores.province);
prov = splitapply(@(x)x(1), WaveScores.province, G);

% 日志统计
n_total     = splitapply(@numel, WaveScores.RRI_wave, G);
n_nan_rri   = splitapply(@(x)sum(isnan(x)), WaveScores.RRI_wave, G);
n_nan_mei   = splitapply(@(x)sum(isnan(x)), WaveScores.MEI_wave, G);
sum_w       = splitapply(@(w)sum(w,'omitnan'), WaveScores.w_tail, G);

% 计算聚合分
RRI_bar = splitapply(@(s,w) aggFun(s,w,A.nan_mode,'RRI'), WaveScores.RRI_wave, WaveScores.w_tail, G);
MEI_bar = splitapply(@(s,w) aggFun(s,w,A.nan_mode,'MEI'), WaveScores.MEI_wave, WaveScores.w_tail, G);

% 记录“有效权重和”
sum_w_used_rri = splitapply(@(s,w) sum(w(~isnan(s))), WaveScores.RRI_wave, WaveScores.w_tail, G);
sum_w_used_mei = splitapply(@(s,w) sum(w(~isnan(s))), WaveScores.MEI_wave, WaveScores.w_tail, G);

ProvAgg = table(prov, n_total, n_nan_rri, n_nan_mei, sum_w, sum_w_used_rri, sum_w_used_mei, ...
                RRI_bar, MEI_bar, ...
                'VariableNames', {'province','n_waves','n_nan_RRI','n_nan_MEI','sum_w_tail','sum_w_used_RRI','sum_w_used_MEI', ...
                                  'RRI_bar','MEI_bar'});

% 排名
ProvAgg.RRI_bar_rank = tiedrank(-ProvAgg.RRI_bar);
ProvAgg.MEI_bar_rank = tiedrank(-ProvAgg.MEI_bar);

% ---- 落盘 ----
outdir = char(A.outdir); if ~isfolder(outdir), mkdir(outdir); end
tag = datestr(now,'yyyymmdd_HHMM');

if A.timestamp
    suffix = "_" + string(tag);
else
    suffix = "";
end

mfile = fullfile(outdir, sprintf('province_agg%s.mat', suffix));
cfile = fullfile(outdir, sprintf('province_agg%s.csv', suffix));

save(mfile, 'ProvAgg');
writetable(ProvAgg, cfile);

if A.verbose
    fprintf('[aggregate_province_scores] Saved: %s\n', cfile);
    fprintf('  nan_mode=%s | 省份数=%d | 中位n_waves=%.1f\n', A.nan_mode, height(ProvAgg), median(ProvAgg.n_waves));
end
end

% ===== helpers =====
function out = aggFun(score, w, nan_mode, tag)
m = ~isnan(score);
switch lower(string(nan_mode))
    case 'exclude'      % 主分析：剔除缺测，权重重归一
        if ~any(m), out = NaN; return; end
        ww = w(m); ww = ww / sum(ww);
        out = sum( ww .* score(m) );
    case 'prov_mean'    % 敏感性：省内均值填补
        mu = mean(score(m),'omitnan');
        score(~m) = mu;
        ww = w / sum(w);
        out = sum( ww .* score );
    otherwise
        error('Unknown nan_mode: %s', nan_mode);
end
end
