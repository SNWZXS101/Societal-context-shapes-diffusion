function ProvProgress = compute_lpi(WaveScores, ProvScores, varargin)
% COMPUTE_LPI
% 基于波级 SRS_wave 的首末差 + 趋势斜率，计算进步指数 LPI（并输出可视化友好表）。
%
% 用法：
%   ProvProgress = compute_lpi(WaveScores, ProvScores, ...
%       'features_mat','features_aug.mat', ...
%       'delta_ceiling',true, ...
%       'slope_method','theilsen', ...    % 'theilsen' | 'robustfit'
%       'w_delta',0.5, 'w_slope',0.5, ...
%       'outdir','scores','timestamp',true,'verbose',true);

p = inputParser; p.FunctionName = mfilename;
addParameter(p,'features_mat','',@(s)ischar(s)||isstring(s));
addParameter(p,'delta_ceiling',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'slope_method','theilsen',@(s)any(strcmpi(s,{'theilsen','robustfit'})));
addParameter(p,'w_delta',0.5,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'w_slope',0.5,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'outdir','scores',@(x)ischar(x)||isstring(x));
addParameter(p,'timestamp',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'verbose',true,@(x)islogical(x)||ismember(x,[0 1]));
parse(p,varargin{:});
A = p.Results;

W = WaveScores;

% ---- 需要列检查 ----
need = {'province','wave_id','variant','RRI_wave','MEI_wave'};
for i = 1:numel(need)
    if ~ismember(need{i}, W.Properties.VariableNames)
        error('WaveScores 缺少列：%s', need{i});
    end
end

% ---- 计算波级综合分 SRS_wave ----
W.SRS_wave = 0.4*W.RRI_wave + 0.6*W.MEI_wave;

% ---- 组装波序（优先 wave_idx → start_date → wave_id 数值化 → 行序）----
order_col = 'row_order';
W.row_order = (1:height(W)).';
if ~isstring(W.province), W.province = string(W.province); end
if ~isstring(W.wave_id),  W.wave_id  = string(W.wave_id);  end
if ~isstring(W.variant),  W.variant  = string(W.variant);  end

% 尝试从 features_aug 取 wave_idx / start_date
if ~isempty(A.features_mat) && isfile(A.features_mat)
    S = load(A.features_mat);
    Tfeat = [];
    cand = {'Features','Tfeat','features','T'};
    for nm = cand
        if isfield(S, nm{1}), Tfeat = S.(nm{1}); break; end
    end
    if istable(Tfeat)
        Tkey = Tfeat(:, intersect({'province','wave_id','wave_idx','start_date'}, Tfeat.Properties.VariableNames));
        if ~isstring(Tkey.province), Tkey.province = string(Tkey.province); end
        if ~isstring(Tkey.wave_id),  Tkey.wave_id  = string(Tkey.wave_id);  end
        W = outerjoin(W, Tkey, 'Keys',{'province','wave_id'}, 'MergeKeys',true);
        if ismember('wave_idx', W.Properties.VariableNames)
            W.order_wave = double(W.wave_idx); order_col='order_wave';
        elseif ismember('start_date', W.Properties.VariableNames)
            try
                W.order_wave = datenum(W.start_date); order_col='order_wave';
            catch
                % 如果 start_date 不是 datetime，尝试转换
                W.order_wave = str2double(string(W.start_date));
                order_col='order_wave';
            end
        end
    end
end

% 若还没有 order_wave，尝试 wave_id 数值化
if ~ismember('order_wave', W.Properties.VariableNames) || all(isnan(W.order_wave))
    tmp = str2double(W.wave_id);
    if any(~isnan(tmp))
        W.order_wave = tmp;
        order_col = 'order_wave';
    else
        % 退路：保留原行序（组内排序时用 row_order）
        W.order_wave = W.row_order;
        order_col = 'order_wave';
    end
end

% ---- 生成省内轨迹（排序）并输出可视化表 ----
Wtraj = sortrows(W(:, {'province','wave_id','variant','SRS_wave','RRI_wave','MEI_wave','order_wave'}), {'province','order_wave'});
Wtraj.wave_seq = within_group_seq(Wtraj.province);   % 1..n（省内）
ProvTrajectories = Wtraj(:, {'province','wave_seq','wave_id','variant','SRS_wave','RRI_wave','MEI_wave'});
ProvTrajectories.Properties.VariableNames{'wave_seq'} = 'order';

% ---- 逐省计算首末差 + 趋势斜率 ----
G = findgroups(ProvTrajectories.province);
prov = splitapply(@(x)x(1), ProvTrajectories.province, G);

first_SRS = splitapply(@(y) pick_first(y), ProvTrajectories.SRS_wave, G);
last_SRS  = splitapply(@(y) pick_last(y),  ProvTrajectories.SRS_wave, G);
delta_raw = last_SRS - first_SRS;

% 天花板/地板校正（将正向/负向改进按可用空间归一）
if A.delta_ceiling
    delta_norm = splitapply(@(y) delta_ceiling_norm(y), ProvTrajectories.SRS_wave, G);
else
    delta_norm = delta_raw;
end

% 趋势斜率
slope = splitapply(@(ord,y) slope_fun(ord,y,A.slope_method), ProvTrajectories.order, ProvTrajectories.SRS_wave, G);

% 标准化到 0–100 并合成 LPI
comp_delta = zto100(delta_norm);
comp_slope = zto100(slope);
LPI = A.w_delta*comp_delta + A.w_slope*comp_slope;

ProvProgress = table(prov, first_SRS, last_SRS, delta_raw, delta_norm, slope, comp_delta, comp_slope, LPI, ...
                     'VariableNames', {'province','first_SRS','last_SRS','delta_SRS','delta_SRS_norm','slope','comp_delta','comp_slope','LPI'});
ProvProgress.LPI_rank = tiedrank(-ProvProgress.LPI);

% ---- 落盘 ----
outdir = char(A.outdir); if ~isfolder(outdir), mkdir(outdir); end
tag = datestr(now,'yyyymmdd_HHMM');

if A.timestamp
    suffix = "_" + string(tag);
else
    suffix = "";
end

mfile = fullfile(outdir, sprintf('province_progress%s.mat', suffix));
cfile = fullfile(outdir, sprintf('province_progress%s.csv', suffix));
save(mfile, 'ProvProgress');
writetable(ProvProgress, cfile);

% 省内轨迹（用于折线/小提琴/变化带）
writetable(ProvTrajectories, fullfile(outdir, sprintf('province_trajectories%s.csv', suffix)));

% 地图友好总表（拼 SRS 与 LPI）
MapReady = innerjoin(ProvScores(:, {'province','SRS'}), ProvProgress(:, {'province','LPI'}));
writetable(MapReady, fullfile(outdir, sprintf('map_ready_SRS_LPI%s.csv', suffix)));

% 排名 Top10 辅助表
Top10_SRS = sortrows(ProvScores(:, {'province','SRS'}),'SRS','descend'); Top10_SRS = Top10_SRS(1:min(10,height(Top10_SRS)),:);
Top10_LPI = sortrows(ProvProgress(:, {'province','LPI'}),'LPI','descend'); Top10_LPI = Top10_LPI(1:min(10,height(Top10_LPI)),:);
writetable(Top10_SRS, fullfile(outdir, sprintf('top10_SRS%s.csv', suffix)));
writetable(Top10_LPI, fullfile(outdir, sprintf('top10_LPI%s.csv', suffix)));

if A.verbose
    fprintf('[compute_lpi] Saved: %s\n', cfile);
    fprintf('  slope=%s | w_delta=%.2f w_slope=%.2f | delta_ceiling=%d\n', A.slope_method, A.w_delta, A.w_slope, A.delta_ceiling);
end

end




% ===== helpers =====
function y = pick_first(v)
v = v(~isnan(v));
if isempty(v), y = NaN; else, y = v(1); end
end
function y = pick_last(v)
v = v(~isnan(v));
if isempty(v), y = NaN; else, y = v(end); end
end
function seq = within_group_seq(g)
seq = zeros(size(g));
u = unique(g,'stable');
for i = 1:numel(u)
    idx = find(g==u(i));
    seq(idx) = 1:numel(idx);
end
end
function d = delta_ceiling_norm(y)
% y: 省内顺序的 SRS_wave 序列
y = y(~isnan(y));
if numel(y)<2, d = NaN; return; end
first = y(1); last = y(end); delta = last - first;
if delta>=0
    cap = max(1, 100 - first);
    d = 100 * delta / cap;    % 0..+100（当 last 达到 100 时为 +100）
else
    cap = max(1, first);
    d = 100 * delta / cap;    % 0..-100（当 last 达到 0 时为 -100）
end
end
function s = slope_fun(ord, y, method)
m = ~isnan(ord) & ~isnan(y);
ord = ord(m); y = y(m);
if numel(y) < 2, s = NaN; return; end
switch lower(string(method))
    case 'theilsen'
        s = theilsen(ord, y);
    case 'robustfit'
        b = robustfit(ord, y); % y = b(1) + b(2)*x
        s = b(2);
    otherwise
        s = NaN;
end
end
function s = theilsen(x,y)
% Theil–Sen 中位数斜率
x = x(:); y = y(:);
n = numel(x);
S = [];
for i = 1:n-1
    dx = x(i+1:end) - x(i);
    dy = y(i+1:end) - y(i);
    m  = dx~=0;
    if any(m)
        S = [S; dy(m)./dx(m)]; %#ok<AGROW>
    end
end
if isempty(S), s = NaN; else, s = median(S,'omitnan'); end
end
function y = zto100(x)
x = double(x(:));
if all(isnan(x))
    y = nan(size(x));
    return;
end
xmin = min(x,[],'omitnan'); xmax = max(x,[],'omitnan');
if ~isfinite(xmin) || ~isfinite(xmax) || abs(xmax-xmin)<1e-12
    y = 50 + 0*x;
else
    y = 100 * (x - xmin) / (xmax - xmin);
end
end
