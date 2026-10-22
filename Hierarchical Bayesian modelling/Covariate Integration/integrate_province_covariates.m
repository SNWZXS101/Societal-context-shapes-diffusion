function integrate_province_covariates(varargin)
% INTEGRATE_PROVINCE_COVARIATES
% 将省级协变量（Excel 30行）合并进 features.mat 的 210 波记录，并输出 features_aug.mat。
% 可选：从 smf_fpca_cluster_results.mat 并入 FPCA xi1/xi2 与 cluster4。
%
% 用法示例：
% integrate_province_covariates( ...
%   'features_mat','features.mat', ...
%   'province_xlsx','province_covariates.xlsx', ...
%   'fpca_mat','smf_fpca_cluster_results.mat', ...
%   'out_mat','features_aug.mat', ...
%   'write_preview_csv', true);

% ----------------------- Parse args -----------------------
p = inputParser;
addParameter(p,'features_mat','features.mat', @ischar);
addParameter(p,'province_xlsx','province_covariates.xlsx', @ischar);
addParameter(p,'fpca_mat','', @ischar);           % 可留空
addParameter(p,'out_mat','features_aug.mat', @ischar);
addParameter(p,'write_preview_csv', false, @islogical);
parse(p,varargin{:});
args = p.Results;

% ----------------------- Load features.mat ----------------
if ~isfile(args.features_mat)
    error('找不到 features_mat: %s', args.features_mat);
end
S = load(args.features_mat);

% 兼容不同变量名（首选 Features）
if isfield(S,'Features')
    T = S.Features;
elseif isfield(S,'Tfeat')
    T = S.Tfeat;
else
    % 尝试寻找 table 变量
    fn = fieldnames(S);
    idx = find(structfun(@istable,S),1);
    if isempty(idx), error('features.mat 中未找到 table 类型变量。'); end
    T = S.(fn{idx});
end

if ~istable(T)
    error('载入的特征数据不是 table。');
end

% 关键列检查
req_keys = {'province','wave_id'};
for k = 1:numel(req_keys)
    if ~ismember(req_keys{k}, T.Properties.VariableNames)
        error('特征表缺少关键列: %s', req_keys{k});
    end
end

% ----------------------- Read province_xlsx ---------------
if ~isfile(args.province_xlsx)
    error('找不到省级协变量文件: %s', args.province_xlsx);
end

try
    provTbl = readtable(args.province_xlsx,'Sheet','province_covars','PreserveVariableNames',true);
catch
    % 兼容没有指定 sheet 的情况
    provTbl = readtable(args.province_xlsx,'PreserveVariableNames',true);
end

if height(provTbl)==0
    error('province_covariates.xlsx 的表为空。请确认文件内容。');
end
if ~ismember('province', provTbl.Properties.VariableNames)
    error('province_covars 必须包含列 "province"。');
end

% ----------------------- Drop possible units row ----------
% 若最后一行是“单位/units”备注，自动删除
lastRow = provTbl(end,:);
isUnitsRow = false;
try
    provVal = "";
    if ismember('province', provTbl.Properties.VariableNames)
        provVal = strtrim(string(lastRow.province));
    end
    % 判断标准：最后一行大量非空文本，且province列为空或包含“单位/units”
    nonEmptyCount = sum(varfun(@(s) sum(~ismissing(s)), lastRow, 'OutputFormat','uniform')>0);
    if nonEmptyCount >= 3
        pv = lower(string(provVal));
        if any(contains(pv, ["单位","unit"]))
            isUnitsRow = true;
        end
    end
catch
    % 忽略异常
end
if isUnitsRow
    provTbl(end,:) = [];
end

% ----------------------- Normalize keys -------------------
T.province       = strtrim(string(T.province));
provTbl.province = strtrim(string(provTbl.province));

% 检查覆盖情况
uProv = unique(T.province);
missingProv = setdiff(uProv, unique(provTbl.province));
if ~isempty(missingProv)
    warning('以下省份在省级协变量中缺失，将产生 NaN：\n  %s', strjoin(missingProv, ', '));
end

% ----------------------- Broadcast join -------------------
% 将省级列（除province）广播到每一波
joinVars = setdiff(provTbl.Properties.VariableNames, {'province'});
% 先删除主表中同名列，避免冲突
keepVars = T.Properties.VariableNames;
overlap = intersect(joinVars, keepVars);
if ~isempty(overlap)
    T(:,overlap) = [];
end
% 为所有新列预创建变量
for j = 1:numel(joinVars)
    v = joinVars{j};
    col = provTbl.(v);
    if iscell(col) || isstring(col) || iscategorical(col)
        T.(v) = strings(height(T),1);
    else
        T.(v) = nan(height(T),1);
    end
end

% 通过映射快速填充
[~,loc] = ismember(T.province, provTbl.province);
for j = 1:numel(joinVars)
    v = joinVars{j};
    src = provTbl.(v);
    if iscell(src) || isstring(src) || iscategorical(src)
        outcol = strings(height(T),1);
        idx = loc>0 & ~ismissing(src(loc));
        outcol(idx) = string(src(loc(idx)));
        T.(v) = outcol;
    else
        outcol = nan(height(T),1);
        idx = loc>0;
        outcol(idx) = double(src(loc(idx)));
        T.(v) = outcol;
    end
end

% ----------------------- Optional: merge FPCA -------------
xi_added = false;
if ~isempty(args.fpca_mat) && isfile(args.fpca_mat)
    F = load(args.fpca_mat);
    scores = [];
    cluster4 = [];
    if isfield(F,'scores'), scores = F.scores; end
    if isfield(F,'cluster')
        if isstruct(F.cluster) && isfield(F.cluster,'labels')
            cluster4 = F.cluster.labels;
        else
            cluster4 = F.cluster;
        end
    elseif isfield(F,'labels')
        cluster4 = F.labels;
    end

    if ~ismember('wave_idx', T.Properties.VariableNames)
        warning('特征表不含 wave_idx，无法按索引并入 xi1/xi2/cluster4。');
    else
        wi = double(T.wave_idx);
        if ~isempty(scores) && size(scores,1)>=2
            T.xi1 = nan(height(T),1);
            T.xi2 = nan(height(T),1);
            good = wi>=1 & wi<=size(scores,2);
            T.xi1(good) = scores(1, wi(good));
            T.xi2(good) = scores(2, wi(good));
            xi_added = true;
        end
        if ~isempty(cluster4)
            T.cluster4 = nan(height(T),1);
            good = wi>=1 & wi<=numel(cluster4);
            T.cluster4(good) = cluster4(wi(good));
        end
    end
end

% ----------------------- Basic transforms ----------------
% 为下游建模准备对数列（不覆盖原列）
numCols = intersect( ...
   {'gdp_pc_2019','pop_density_2020','total_density_2020', ...
    'hospital_beds_per1k_2019','passenger_volume_2019', ...
    'transit_volume_2019','inbound_visitors_2019', ...
    'inc_disp_pc_2019','inc_disp_pc_urban_2019','inc_disp_pc_rural_2019'}, ...
    T.Properties.VariableNames);

for j = 1:numel(numCols)
    v = numCols{j};
    lv = "log_"+v;
    if ~ismember(lv, T.Properties.VariableNames)
        x = double(T.(v));
        x(~isfinite(x) | x<=0) = NaN;   % 避免非正数取对数
        T.(lv) = log(x);
    end
end

% ----------------------- Save & report --------------------
Features = T; %#ok<NASGU>
Meta = struct();
Meta.generated_on = datestr(now,'yyyy-mm-dd HH:MM:SS');
Meta.source_features_mat = args.features_mat;
Meta.province_xlsx = args.province_xlsx;
Meta.fpca_mat = args.fpca_mat;
Meta.added_fpca = xi_added;

save(args.out_mat,'Features','Meta');
fprintf('[OK] 已保存合并结果到 %s. 行=%d, 列=%d\n', args.out_mat, height(T), width(T));

if args.write_preview_csv
    try
        outcsv = strrep(args.out_mat,'.mat','_preview.csv');
        writetable(T, outcsv);
        fprintf('已输出预览 CSV: %s\n', outcsv);
    catch ME
        warning('写预览CSV失败：%s', ME.message);
    end
end
end
