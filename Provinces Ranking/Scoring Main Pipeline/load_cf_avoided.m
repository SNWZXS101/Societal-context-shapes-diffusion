function Tm = load_cf_avoided(cf_results_path, scenario_name)
% LOAD_CF_AVOIDED  读取并筛选指定场景的 avoided_pct，输出两列主键 + pct_avoided_main
% cf_results_path: 指向 cf_results.mat 或 cf_results.csv 或其所在目录
% scenario_name:   例如 'TGC_R2ZERO' / 'TGC_STRONG' / ...

if isfolder(cf_results_path)
    matf = fullfile(cf_results_path, 'cf_results.mat');
    csvf = fullfile(cf_results_path, 'cf_results.csv');
else
    [p,f,e] = fileparts(cf_results_path);
    if strcmpi(e,'.mat'), matf = cf_results_path; csvf = ''; else, matf = ''; csvf = cf_results_path; end
end

Results = [];
if ~isempty(matf) && isfile(matf)
    S = load(matf);
    % 兼容字段名
    if isfield(S,'Results')
        Results = S.Results;
    elseif isfield(S,'ResTbl')
        Results = S.ResTbl;
    end
elseif ~isempty(csvf) && isfile(csvf)
    Results = readtable(csvf);
else
    error('未找到 cf_results 文件：%s', cf_results_path);
end
if ~istable(Results)
    error('cf_results 无法解析为 table');
end

% 统一字符串
for c = ["province","wave_id","scenario"]
    if ismember(c, Results.Properties.VariableNames)
        if ~isstring(Results.(c)), Results.(c) = string(Results.(c)); end
    else
        error('cf_results 缺少列：%s', c);
    end
end
if ~ismember("avoided_pct", Results.Properties.VariableNames)
    error('cf_results 缺少列：avoided_pct');
end

% 筛选场景
mask = upper(strtrim(Results.scenario)) == upper(strtrim(string(scenario_name)));
R = Results(mask, {'province','wave_id','avoided_pct'});
R.Properties.VariableNames{'avoided_pct'} = 'pct_avoided_main';

% 去重（若有重复，取均值）
[G, prov, wid] = findgroups(R.province, R.wave_id);
val = splitapply(@mean, R.pct_avoided_main, G);
Tm = table(prov, wid, val, 'VariableNames', {'province','wave_id','pct_avoided_main'});

end
