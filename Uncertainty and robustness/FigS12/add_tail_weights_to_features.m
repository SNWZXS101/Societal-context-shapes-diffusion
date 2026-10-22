function add_tail_weights_to_features(features_mat, varargin)
% ADD_TAIL_WEIGHTS_TO_FEATURES
% 依据公式：
%   T_tail = T - max(t_fmax, t_inc_slope_max), 截断为 >=0
%   w_tail_k = min(1, T_tail / k), k ∈ {7,10,14}
% 在 features_mat 的表 Features（或 Tfeat/Features 等同名备选）中
% 增加/更新列：w_tail_7, w_tail_10, w_tail_14。
%
% 用法：
%   add_tail_weights_to_features('features_aug.mat');                 % 覆盖保存
%   add_tail_weights_to_features('features_aug.mat','save_as','features_aug_tail.mat'); % 另存为
%
% 要求 Features 表至少包含列：T, t_fmax, t_inc_slope_max（区分大小写不严）。

p = inputParser; p.FunctionName = mfilename;
addRequired(p,'features_mat',@(s)ischar(s)||isstring(s));
addParameter(p,'save_as','',@(s)ischar(s)||isstring(s));  % 为空则覆盖原文件
parse(p,features_mat,varargin{:});
A = p.Results;

if ~isfile(A.features_mat)
    error('找不到文件：%s', A.features_mat);
end
S = load(A.features_mat);
Tfeat = [];
cands = {'Features','Tfeat','features','T'};
for k = 1:numel(cands)
    if isfield(S, cands{k})
        Tfeat = S.(cands{k});
        feat_name = cands{k};
        break;
    end
end
if ~istable(Tfeat)
    error('features_mat 中未找到表（Features/Tfeat/features/T）');
end

% 兼容大小写和下划线差异：做一层宽松匹配
function col = pick(cols, names)
    col = '';
    for i=1:numel(names)
        m = strcmpi(cols, names{i});
        if any(m), col = cols{find(m,1)}; return; end
    end
end

cols = Tfeat.Properties.VariableNames;
c_T   = pick(cols, {'T','t_len','tL','t_l','t_total'});
c_tf  = pick(cols, {'t_fmax','tf2_peak','t_tf2_peak','t_f_max'});
c_inc = pick(cols, {'t_inc_slope_max','t_inc_slope','t_incmax','t_inc_max'});

need_msg = '';
if isempty(c_T),   need_msg = [need_msg, ' T']; end
if isempty(c_tf),  need_msg = [need_msg, ' t_fmax']; end
if isempty(c_inc), need_msg = [need_msg, ' t_inc_slope_max']; end
if ~isempty(need_msg)
    error('缺少必要列：%s（请确认列名）', strtrim(need_msg));
end

% 计算 T_tail
T  = double(Tfeat.(c_T));
tf = double(Tfeat.(c_tf));
ti = double(Tfeat.(c_inc));
m  = max(tf, ti);
T_tail = T - m;
T_tail(~isfinite(T_tail)) = NaN;
T_tail(T_tail < 0) = 0;

% 计算 w_tail_k
wt7  = min(1, T_tail/7);
wt10 = min(1, T_tail/10);
wt14 = min(1, T_tail/14);

% 写入/更新列
Tfeat.w_tail_7  = wt7;
Tfeat.w_tail_10 = wt10;
Tfeat.w_tail_14 = wt14;

% 回存
S.(feat_name) = Tfeat;
save_path = char(A.features_mat);
if ~isempty(A.save_as)
    save_path = char(A.save_as);
end
save(save_path, '-struct', 'S');

fprintf('[add_tail_weights_to_features] 已写入列 w_tail_7/10/14 → %s（表名：%s）\n', save_path, feat_name);
% 小提示：可在 compute_wave_scores 中让 w_tail 指向哪个窗口；或在 UQ 中用三窗口做稳健性。
end
