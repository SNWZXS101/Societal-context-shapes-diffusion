function prep_hbm_inputs(varargin)
% PREP_HBM_INPUTS (parametrized + row dropping + abslog + mean-impute)
% 读取 features_aug.mat（及可选 FPCA），完成：
%   - xi1/xi2 灌入（可选符号对齐）
%   - 连续协变量：可选对数 → 标准化 → 均值插补(把 z 的 NaN 置 0)
%   - 类别 one-hot（region/climate；参考 East/Temperate）
%   - 目标变量：identity / log / asinh / abslog（|x|取log）
%   - 对指定目标的非正值行“硬删除”（整行剔除）
%   - 剔除全 NaN 的连续列
%   - 缓存 + 预处理口径保存（记录变换、插补与行删）

% ===== 参数 =====
p = inputParser;
addParameter(p,'features_mat','features_aug.mat',@ischar);
addParameter(p,'fpca_mat','smf_fpca_cluster_results.mat',@ischar);
addParameter(p,'out_dir','results_hbm',@ischar);

% 目标/协变量筛选
addParameter(p,'keep_targets',{},@(c)iscellstr(c)||isstring(c));
addParameter(p,'drop_targets',{},@(c)iscellstr(c)||isstring(c));
addParameter(p,'keep_covars',{},@(c)iscellstr(c)||isstring(c));
addParameter(p,'drop_covars',{},@(c)iscellstr(c)||isstring(c));

% FPC 符号对齐
addParameter(p,'align_xi',false,@islogical);
addParameter(p,'xi_anchor','delta_early',@ischar);

% 目标变换
addParameter(p,'log_targets',{},@(c)iscellstr(c)||isstring(c));
addParameter(p,'asinh_targets',{},@(c)iscellstr(c)||isstring(c));
addParameter(p,'abslog_targets',{},@(c)iscellstr(c)||isstring(c));
addParameter(p,'identity_targets',{},@(c)iscellstr(c)||isstring(c));
addParameter(p,'hard_drop_nonpos_targets',{},@(c)iscellstr(c)||isstring(c));
addParameter(p,'log_nonpos_to_nan',true,@islogical);
addParameter(p,'log_warn_on_nonpos',true,@islogical);

% 连续协变量插补（新增开关，默认开）
addParameter(p,'impute_cont_to_mean',true,@islogical);

parse(p,varargin{:});
args = p.Results;

if ~isfolder(args.out_dir), mkdir(args.out_dir); end
cache_dir = fullfile(args.out_dir,'cache');
if ~isfolder(cache_dir), mkdir(cache_dir); end

% ===== 加载主表 =====
S = load(args.features_mat);
if isfield(S,'Features'), T = S.Features; else, error('features_aug.mat 中未找到 table 变量 "Features"。'); end
if ~istable(T), error('Features 不是 table'); end

% ===== 关键列检查 =====
need_cols = {'province','variant','wave_idx'};
for k=1:numel(need_cols)
    if ~ismember(need_cols{k}, T.Properties.VariableNames)
        error('Features 缺少关键列：%s', need_cols{k});
    end
end
T.province = string(strtrim(T.province));
T.variant  = string(strtrim(T.variant));

% ===== 协变量清单（允许缺列补 NaN） =====
covar_orig_all = { ...
 'age65 share', ...
 'climate zone', ...
 'education share', ...
 'gdp_{pc}', ...
 'hospital beds per 1k', ...
 'proportion of students', ...  
 'inbound overnight visitors', ...
 'disposable income_{pc}', ...
 'passenger volume', ...
 'public transit volume', ...
 'region category', ...
 'tertiary-sector share', ...
 'urban population density', ...
 'urban rate' };

for i=1:numel(covar_orig_all)
    if ~ismember(covar_orig_all{i}, T.Properties.VariableNames)
        warning('协变量缺失：%s（将用 NaN）', covar_orig_all{i});
        T.(covar_orig_all{i}) = NaN(height(T),1);
    end
end

% ===== 目标集合（keep_targets 覆盖默认） =====
targets_default = {'delta_early','kappa_cross','rho_eps','tf1_mean_early','tf2_peak','t_tf2_peak'};
if isempty(args.keep_targets), targets = targets_default;
else, targets = cellstr(string(args.keep_targets));
end
if ~isempty(args.drop_targets)
    targets = setdiff(targets, cellstr(string(args.drop_targets)), 'stable');
end
for i=1:numel(targets)
    if ~ismember(targets{i}, T.Properties.VariableNames)
        T.(targets{i}) = NaN(height(T),1);
    end
end

% ===== 行剔除：对指定目标的非正值整行删除 =====
drop_set = lower(string(args.hard_drop_nonpos_targets));
drop_mask = false(height(T),1);
for i=1:numel(drop_set)
    nm = char(drop_set(i));
    if ismember(nm, lower(string(T.Properties.VariableNames)))
        y = double(T.(nm));
        bad = isfinite(y) & (y<=0);
        drop_mask = drop_mask | bad;
    end
end
n_drop = nnz(drop_mask);
if n_drop>0
    T = T(~drop_mask, :);
    fprintf('[prep] 硬删除行：%d 条（原因：%s 的非正值）\n', n_drop, strjoin(cellstr(drop_set), ', '));
end

% ===== 从 FPCA 结果按 wave_idx 灌入 xi1/xi2（若存在） =====
fpca_ok = (exist(args.fpca_mat,'file')==2);
if fpca_ok
    try
        Fp = load(args.fpca_mat);
        Sscores = [];
        if isstruct(Fp) && isfield(Fp,'scores')
            Sscores = Fp.scores;
        elseif isstruct(Fp)
            fn = fieldnames(Fp);
            for k=1:numel(fn)
                if isstruct(Fp.(fn{k})) && isfield(Fp.(fn{k}),'scores')
                    tmp = Fp.(fn{k}).scores;
                    if ~isempty(tmp), Sscores = tmp; break; end
                end
            end
        end

        if ~isempty(Sscores)
            Nnow = height(T);
            [r,c] = size(Sscores);
            if ~ismember('xi1', T.Properties.VariableNames), T.xi1 = NaN(Nnow,1); end
            if ~ismember('xi2', T.Properties.VariableNames), T.xi2 = NaN(Nnow,1); end

            wi = double(T.wave_idx);
            if r==Nnow && c>=2
                inrng = isfinite(wi) & wi>=1 & wi<=r;
                T.xi1(inrng) = Sscores(wi(inrng),1);
                T.xi2(inrng) = Sscores(wi(inrng),2);
                left = isnan(T.xi1) | isnan(T.xi2);
                if any(left)
                    idx = find(left); idx = idx(idx<=r);
                    T.xi1(idx) = Sscores(idx,1);
                    T.xi2(idx) = Sscores(idx,2);
                end
            elseif c==Nnow && r>=2
                inrng = isfinite(wi) & wi>=1 & wi<=c;
                T.xi1(inrng) = Sscores(1, wi(inrng));
                T.xi2(inrng) = Sscores(2, wi(inrng));
                left = isnan(T.xi1) | isnan(T.xi2);
                if any(left)
                    idx = find(left); idx = idx(idx<=c);
                    T.xi1(idx) = Sscores(1, idx);
                    T.xi2(idx) = Sscores(2, idx);
                end
            else
                warning('prep_hbm_inputs:fpcaSizeMismatch', ...
                    'FPCA scores 尺寸为 (%d,%d) 与当前 N=%d 不匹配，跳过自动映射。', r, c, Nnow);
            end

            if args.align_xi && ismember(args.xi_anchor, T.Properties.VariableNames)
                anchor = double(T.(args.xi_anchor));
                s1 = sign(nansum(T.xi1 .* anchor)); if s1==0, s1=1; end
                s2 = sign(nansum(T.xi2 .* anchor)); if s2==0, s2=1; end
                T.xi1 = s1*T.xi1; T.xi2 = s2*T.xi2;
            end
        end
    catch err
        warning('prep_hbm_inputs:fpcaLoad', '加载 FPCA 结果失败（已忽略）：%s', err.message);
    end
end

% ===== 名称映射 =====
safe = @(s) regexprep(regexprep(s,'\s+','_'), '[^A-Za-z0-9_]', '');
covar_all = covar_orig_all;
if ~isempty(args.keep_covars)
    covar_all = intersect(string(covar_all), string(args.keep_covars), 'stable');
end
if ~isempty(args.drop_covars)
    covar_all = setdiff(string(covar_all), string(args.drop_covars), 'stable');
end
covar_all = cellstr(covar_all);

covar_safe = cellfun(safe, covar_all, 'UniformOutput', false);
is_factor = ismember(covar_all, {'region category','climate zone'});
cont_orig = covar_all(~is_factor);
cont_safe = covar_safe(~is_factor);
fac_orig  = covar_all(is_factor);

% ===== 连续协变量是否先 log（用于“先 log 再标准化”） =====
normkey = @(s) lower(regexprep(s,'[^a-zA-Z0-9]',''));
is_log_col = false(1,numel(cont_orig));
for j=1:numel(cont_orig)
    nm = normkey(cont_orig{j});
    if (contains(nm,'gdp') && contains(nm,'pc')) || ...
       (contains(nm,'income') && contains(nm,'pc')) || ...
       strcmp(nm,'urbanpopulationdensity') || ...
       contains(nm,'passengervolume') || ...
       contains(nm,'publictransitvolume') || ...
       contains(nm,'inboundovernightvisitors') || ...
       contains(nm,'hospitalbedsper1k')
        is_log_col(j) = true;
    end
end

% ===== 连续协变量：变换 + 标准化 + 均值插补 =====
stats = struct();
X_cont = nan(height(T), numel(cont_orig));
miss_rate = zeros(1,numel(cont_orig));
for j=1:numel(cont_orig)
    c = cont_orig{j};
    x = double(T.(c));
    if is_log_col(j)
        x(~isfinite(x) | x<=0) = NaN; x = log(x);
        stats.(safe(['log_' c])).transform = 'log';
    else
        stats.(safe(c)).transform = 'identity';
    end
    mu = nanmean(x); sd = nanstd(x);
    if sd==0 || ~isfinite(sd), sd=1; end
    z = (x - mu) ./ sd;

    % —— 均值插补：把 z 中的 NaN 置 0（避免后续行被全剔除）
    miss = ~isfinite(z);
    miss_rate(j) = mean(miss);
    if args.impute_cont_to_mean
        z(miss) = 0;
    end

    X_cont(:,j) = z;
    stats.(safe(c)).mean = mu;
    stats.(safe(c)).sd   = sd;
    stats.(safe(c)).imputed_to_zero = args.impute_cont_to_mean;
    stats.(safe(c)).missing_rate = miss_rate(j);
end

% —— 剔除“全为缺失（插补前全 NaN）”的连续列 —— 
allnan_cont = miss_rate == 1;
if any(allnan_cont)
    bad_names = cont_orig(allnan_cont);
    warning('剔除全 NaN 的连续协变量列：%s', strjoin(bad_names, ', '));
    X_cont(:, allnan_cont) = [];
    cont_orig(allnan_cont) = [];
    cont_safe(allnan_cont) = [];
end

% ===== 类别 one-hot（参考 East/Temperate） =====
reg_levels = ["North","Northeast","East","Central","South","Southwest","Northwest"]; % 参考 East
cli_levels = ["Temperate","Subtropical","Tropic","Plateau","Arid"];                   % 参考 Temperate
reg = string(T.("region category")); reg = strtrim(reg);
cli = string(T.("climate zone"));    cli = strtrim(cli);

reg_ref = "East";      reg_cols = setdiff(reg_levels, reg_ref, 'stable');
cli_ref = "Temperate"; cli_cols = setdiff(cli_levels, cli_ref, 'stable');

X_reg = zeros(height(T), numel(reg_cols));
X_cli = zeros(height(T), numel(cli_cols));
for k=1:numel(reg_cols), X_reg(:,k) = double(reg==reg_cols(k)); end
for k=1:numel(cli_cols), X_cli(:,k) = double(cli==cli_cols(k)); end

% ===== 省/变异株索引 =====
[prov_names, ~, prov_id] = unique(T.province, 'stable');
var_levels = ["Original","Delta","Omicron"];
var_map = zeros(height(T),1);
for i=1:height(T)
    vi = T.variant(i); idx = find(var_levels==vi, 1);
    if isempty(idx), idx = 1; end
    var_map(i) = idx;
end

% ===== 固定效应设计矩阵 =====
X_fe = [X_cont, X_reg, X_cli];
fe_names_orig = [cont_orig, cellstr("region "+reg_cols), cellstr("climate "+cli_cols)];

% ===== 目标变量变换集合（默认 + 参数覆盖） =====
log_defaults    = {'delta_early','rho_eps','tf1_mean_early','tf2_peak','t_tf2_peak'};
abslog_defaults = {'kappa_cross'};
asinh_defaults  = {};

to_lower = @(x) lower(string(x));
log_set    = unique([to_lower(log_defaults),    to_lower(args.log_targets)], 'stable');
abslog_set = unique([to_lower(abslog_defaults), to_lower(args.abslog_targets)], 'stable');
asinh_set  = unique([to_lower(asinh_defaults),  to_lower(args.asinh_targets)], 'stable');
ident_set  = to_lower(args.identity_targets);

Y = struct(); tf_info = struct();
for i=1:numel(targets)
    key  = targets{i};
    keyL = lower(string(key));
    y = double(T.(key));

    % 软处理：对 log/abslog 的非正值置 NaN（若未被硬删）
    is_loglike = any(keyL==log_set) || any(keyL==abslog_set);
    if is_loglike
        if any(keyL==log_set)
            badcnt = sum(isfinite(y) & y<=0);
            if badcnt>0 && args.log_warn_on_nonpos
                warning('%s: %d 非正值遇到 log 变换（将置 NaN）。', key, badcnt);
            end
            if args.log_nonpos_to_nan, y(isfinite(y) & y<=0) = NaN; end
        elseif any(keyL==abslog_set)
            badcnt = sum(isfinite(y) & y==0);
            if badcnt>0 && args.log_warn_on_nonpos
                warning('%s: %d 为 0（abslog 变换将置 NaN）。', key, badcnt);
            end
            if args.log_nonpos_to_nan, y(isfinite(y) & y==0) = NaN; end
        end
    end

    % 正式变换（优先级：identity > asinh > abslog > log > identity）
    if any(keyL == ident_set), tf='identity'; y_tr = y;
    elseif any(keyL == asinh_set), tf='asinh'; y_tr = asinh(y);
    elseif any(keyL == abslog_set), tf='abslog'; y_tr = log(abs(y));
    elseif any(keyL == log_set), tf='log'; y_tr = log(y);
    else, tf='identity'; y_tr = y;
    end

    Y.(key) = y_tr;
    tf_info.(key) = tf;
end

% ===== 保存缓存 =====
cache = struct();
cache.X_fe = X_fe;
cache.fe_names_orig = fe_names_orig;
cache.covar_orig = covar_all;
cache.cont_orig = cont_orig;
cache.cont_safe = cont_safe;
cache.reg_cols = cellstr(reg_cols);
cache.cli_cols = cellstr(cli_cols);
cache.prov_id = prov_id;
cache.prov_names = prov_names;
cache.var_id = var_map;
cache.var_levels = cellstr(var_levels);
cache.Y = Y;
cache.tf_info = tf_info;
cache.N = height(T);
cache.K = size(X_fe,2);
cache.P = numel(prov_names);
cache.V = numel(var_levels);
cache.active_fe_names = fe_names_orig;
cache.active_targets  = targets;
cache.row_drop = struct('n_drop', n_drop, 'targets', cellstr(drop_set));
cache.impute_cont_to_mean = args.impute_cont_to_mean;

save(fullfile(cache_dir,'hbm_cache.mat'),'cache');

% ===== 预处理口径（manifest） =====
manifest = struct();
name_map_tbl = table(covar_all(:), cellfun(@(x)regexprep(regexprep(x,'\s+','_'),'[^A-Za-z0-9_]',''), covar_all(:), 'uni',0), ...
    'VariableNames', {'original','safe'});
manifest.name_map = name_map_tbl;
manifest.transform_stats = stats;
manifest.region_levels = cellstr(reg_levels); manifest.region_ref = char(reg_ref);
manifest.climate_levels = cellstr(cli_levels); manifest.climate_ref = char(cli_ref);
manifest.targets = targets;
manifest.target_tf = tf_info;
manifest.row_drop = cache.row_drop;
manifest.impute_cont_to_mean = args.impute_cont_to_mean;
manifest.note = '连续固定效应已标准化；缺失值按均值插补（z=0）；部分连续列先对数；目标变量按 target_tf 变换；hard_drop_nonpos_targets 对应的非正行已整行删除。';

json = jsonencode(manifest,'PrettyPrint',true);
fid = fopen(fullfile(args.out_dir,'preproc_manifest.json'),'w'); fwrite(fid, json); fclose(fid);

% 简短回显 + 每个 target 的可用 N（此时 X_fe 无 NaN）
fprintf('[prep_hbm_inputs+] OK. N=%d, K=%d, P=%d, V=%d | targets=%s\n', ...
    cache.N, cache.K, cache.P, cache.V, strjoin(string(targets),', '));
keys = fieldnames(tf_info); tmp = cell(numel(keys),1);
for i=1:numel(keys), tmp{i} = sprintf('%s->%s', keys{i}, tf_info.(keys{i})); end
fprintf('  target transforms: %s\n', strjoin(tmp, '; '));
for i=1:numel(targets)
    key = targets{i};
    n_ok = sum(isfinite(Y.(key)));
    fprintf('  usable N for %s = %d\n', key, n_ok);
end
end
