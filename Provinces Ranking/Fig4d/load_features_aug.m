function T = load_features_aug(features_mat)
% LOAD_FEATURES_AUG  读取 features_aug.mat 并标准化字段名（鲁棒匹配列名）
% 输出 table T（逐波 210 行），至少包含文中所需关键列

if ~isfile(features_mat)
    error('File not found: %s', features_mat);
end
S = load(features_mat);

% 猜测主 table 名
cand = {'Features','Tfeat','features','T','Features_aug'};
hit = '';
for nm = cand
    if isfield(S, nm{1}), hit = nm{1}; break; end
end
if isempty(hit)
    error('features_aug.mat 未找到主表（候选：%s）', strjoin(cand, ', '));
end
Traw = S.(hit);
if ~istable(Traw)
    error('主对象不是 table：%s', hit);
end
T = Traw;

% ====== 列名鲁棒匹配工具 ======
V = string(T.Properties.VariableNames);
lowerV = lower(V);

    function real = find_col_name(pat_list, fallback_contains)
        % pat_list: cellstr 的正则（lower-case）
        % fallback_contains: 退路匹配（lower-case substring 或空）
        real = "";
        % 1) 正则全量匹配（逐一尝试）
        for ii = 1:numel(pat_list)
            pat = pat_list{ii};
            m = ~cellfun('isempty', regexp(lowerV, pat, 'once'));
            if any(m)
                real = V(find(m,1,'first'));
                return;
            end
        end
        % 2) 退而求其次：contains
        if ~isempty(fallback_contains)
            m = contains(lowerV, fallback_contains);
            if any(m)
                real = V(find(m,1,'first'));
            end
        end
    end

% ====== 关键列名映射（宽松）======
map = struct();

% 主键
map.province = find_col_name({'^province$','^prov(ince)?$'}, 'province');
map.wave_id  = find_col_name({'^wave[_]?id$','^waveid$','^wave$'}, 'wave');

% 权重列（w_tail_14 的各种写法）
map.w_tail = find_col_name({ ...
    '^w[_]?tail[_]?14$', ...          % w_tail_14 / wtail_14 / w_tail14
    '^w[_]?tail.*14$', ...            % w_tailXXXX14
    '^w[_]?tail$' ...                 % 退路：只有 w_tail（不建议，但兼容）
    }, 'w_tail');

% 其他特征
map.delta_early    = find_col_name({'^delta[_]?early$'}, 'delta');
map.kappa_cross    = find_col_name({'^kappa[_]?cross$','^kappa$'}, 'kappa');
map.tf1_mean_early = find_col_name({'^t?f?1[_]?mean[_]?early$','^tf1.*early$'}, 'tf1');
map.t_cross        = find_col_name({'^t[_]?cross$','^tcross$','^t[_]?epsilon$'}, 'cross');
map.rho_eps        = find_col_name({'^rho[_]?eps$','^rho[_]?epsilon$','^rho$'}, 'rho');
map.tf2_peak       = find_col_name({'^t?f?2[_]?peak$','^tf2.*peak$'}, 'tf2');
map.t_tf2_peak     = find_col_name({'^t[_]?tf2[_]?peak$','^ttf2.*peak$'}, 't_tf2');
map.flag_multipart = find_col_name({'^flag[_]?multipart$','^multipart$'}, 'multipart');
map.variant        = find_col_name({'^variant$','^strain$','^lineage$'}, 'variant');

% 应用映射：找到即重命名
keys = fieldnames(map);
for i = 1:numel(keys)
    k = keys{i};
    real = map.(k);
    if strlength(real) > 0 && real ~= k
       % 仅当实际列名存在且与目标不同才重命名
       T.Properties.VariableNames{real} = k;
       fprintf('[load_features_aug] 识别列 %-16s ← %s\n', k, real);
    end
end

% —— 主键存在性检查 ——
req_keys = {'province','wave_id'};
for i = 1:numel(req_keys)
    if ~ismember(req_keys{i}, T.Properties.VariableNames)
        error('缺少关键列：%s', req_keys{i});
    end
end

% —— w_tail 保底/提示 ——
if ~ismember('w_tail', T.Properties.VariableNames)
    warning('未发现 w_tail_14（或其变体），自动回退设为 1');
    T.w_tail = ones(height(T),1);
end

% 统一字符串类型（若存在）
for c = ["province","variant","wave_id"]
    if ismember(c, T.Properties.VariableNames)
        if ~isstring(T.(c)), T.(c) = string(T.(c)); end
    end
end

% flag_multipart 保底
if ~ismember('flag_multipart', T.Properties.VariableNames)
    T.flag_multipart = false(height(T),1);
end

end
