function Features = run_extract_all_features(matPath, varargin)
% RUN_EXTRACT_ALL_FEATURES
% Batch-extract features for all waves in waves_dataset.mat
%
% Usage (examples):
%   Features = run_extract_all_features('waves_dataset.mat');
%   Features = run_extract_all_features('waves_dataset.mat', ...
%       'outMat','features.mat','outCsv','features.csv', ...
%       'useParallel',false,'verbose',true, ...
%       'c_eps',1.2,'q_early',0.30, ...
%       'c_eps_hi',1.0,'c_eps_lo',0.7, ...
%       'alpha_hi',0.15,'alpha_lo',0.10, ...
%       'eps_mode','hybrid-min', ...
%       'm_neg',2,'m_nz',2);

% -------- args --------
p = inputParser; p.FunctionName = mfilename;
addRequired(p,'matPath',@(x)ischar(x)||isstring(x));
addParameter(p,'outMat','features.mat',@(x)ischar(x)||isstring(x));
addParameter(p,'outCsv','features.csv',@(x)ischar(x)||isstring(x));
addParameter(p,'useParallel',false,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'verbose',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'c_eps',0.3,@isnumeric);
addParameter(p,'q_early',0.30,@isnumeric);

% —— 混合阈值与穿越判定的可调参数 —— 
addParameter(p,'c_eps_hi',0.3,@isnumeric);            % MAD上阈系数
addParameter(p,'c_eps_lo',0.3,@isnumeric);            % MAD下阈系数
addParameter(p,'alpha_hi',0.10,@isnumeric);           % f_max上阈比例
addParameter(p,'alpha_lo',0.10,@isnumeric);           % f_max下阈比例
addParameter(p,'eps_mode','hybrid-min',@(x)ischar(x)||isstring(x));  % 'mad'|'fraction'|'hybrid-min'
addParameter(p,'m_neg',1,@(x)isnumeric && isscalar(x) && x>=1);      % 连续负相点
addParameter(p,'m_nz',1,@(x)isnumeric && isscalar(x) && x>=1);       % 连续近零点

parse(p,matPath,varargin{:});
A = p.Results;

% -------- load dataset --------
if ~isfile(matPath), error('File not found: %s', matPath); end
S = load(matPath, 'Waves','Provinces','Meta');
Waves = S.Waves;

% 把 waves_dataset 里的 f 和 init 两列组装为 cell，保持与 Waves 顺序一致
f_cell   = arrayfun(@(w) w.f,   S.Waves, 'UniformOutput', false);
init_cell= arrayfun(@(w) w.init,S.Waves, 'UniformOutput', false);
start_date_vec= arrayfun(@(w) w.start_date, S.Waves); 

N = numel(Waves);
if A.verbose, fprintf('[features] Extracting features for %d waves...\n', N); end

% params pack（传给 compute_features_from_wave）
params = struct( ...
    'c_eps',A.c_eps,'q_early',A.q_early, ...
    'c_eps_hi',A.c_eps_hi,'c_eps_lo',A.c_eps_lo, ...
    'alpha_hi',A.alpha_hi,'alpha_lo',A.alpha_lo, ...
    'eps_mode',string(A.eps_mode), ...
    'm_neg',A.m_neg,'m_nz',A.m_nz );

% prepare outputs
Fcell = cell(N,1);

% Safe parallel switch
usePar = A.useParallel;
try
    v = ver; hasPCT = any(strcmp({v.Name}, 'Parallel Computing Toolbox'));
catch
    hasPCT = false;
end
if usePar && ~hasPCT, usePar = false; end
if usePar
    try
        pool = gcp('nocreate'); if isempty(pool), parpool; end
    catch
        usePar = false;
    end
end

% -------- loop --------
if usePar
    parfor k = 1:N
        Fcell{k} = one_wave_features(Waves(k), params);
    end
else
    for k = 1:N
        Fcell{k} = one_wave_features(Waves(k), params);
    end
end

% Filter successful
ok = ~cellfun(@isempty, Fcell);
Fstruct = [Fcell{ok}]';

% Convert to table
Features = struct2table(Fstruct);

% ===== 覆盖 flag_multipart：以 Excel 为准 =====
try
    % 同目录优先找两个可能的文件名
    cand = {'flag_multiple.xlsx','flag_multipart.xlsx'};
    excelPath = '';
    for ii = 1:numel(cand)
        if isfile(cand{ii})
            excelPath = cand{ii};
            break;
        end
    end

    if ~isempty(excelPath)
        X = readtable(excelPath);   % 需要列：wave_id, flag_multipart（不区分大小写）

        % 宽松匹配列名（显式查索引，避免匿名函数）
        vn = lower(string(X.Properties.VariableNames));
        col_wave = find(ismember(vn, "wave_id"), 1);
        col_flag = find(ismember(vn, "flag_multipart"), 1);
        if isempty(col_wave) || isempty(col_flag)
            error('Excel 必须包含列 wave_id 与 flag_multipart');
        end

        waveX = string(X.(col_wave));
        flagX = double(X.(col_flag));

        % 去重（相同 wave_id 取首个）
        [waveU, ia] = unique(waveX, 'stable');
        flagU = flagX(ia);

        % 左连接覆盖
        Map = table(waveU, flagU, 'VariableNames',{'wave_id','flag_new'});

        if ~ismember('wave_id', Features.Properties.VariableNames)
            error('Features 表缺少 wave_id 列，无法覆盖 flag_multipart');
        end
        Features.wave_id = string(Features.wave_id);

        J = outerjoin(Features, Map, ...
            'Keys','wave_id', 'MergeKeys',true, 'Type','left', ...
            'LeftVariables', Features.Properties.VariableNames, ...
            'RightVariables','flag_new');

        if ~ismember('flag_multipart', Features.Properties.VariableNames)
            Features.flag_multipart = nan(height(Features),1);
        end
        m = ~isnan(J.flag_new);
        Features.flag_multipart(m) = J.flag_new(m);

        if exist('A','var') && isfield(A,'verbose') && A.verbose
            fprintf('[features] Excel override: %d/%d rows (%s)\n', sum(m), height(Features), excelPath);
        end
    else
        if exist('A','var') && isfield(A,'verbose') && A.verbose
            fprintf('[features] Excel override skipped (file not found).\n');
        end
    end
catch ME
    warning(sprintf('Excel 覆盖 flag_multipart 失败：%s', ME.message));
end
% ===== 覆盖结束 =====






% Save
if ~isempty(A.outMat)
    save(A.outMat, 'Features', 'f_cell', 'init_cell', 'start_date_vec', '-v7.3');  % <<< 新增
end


if ~isempty(A.outCsv), writetable(Features, A.outCsv); end

if A.verbose
    fprintf('[features] Done. Rows=%d saved to %s / %s\n', height(Features), A.outMat, A.outCsv);
end
end

% ---------- sub ----------
function F = one_wave_features(W, params)
F = struct();
try
    % ---- meta ----
    F.province  = string(W.province);
    F.wave_idx  = double(W.wave_idx);
    F.wave_id   = string(sprintf('%s-%02d', char(W.province), W.wave_idx));
    F.beta      = double(W.beta);
    F.start_date = W.start_date;

    % ---- 变异株：按 beta 截到 1 位小数判定（0.3/0.7/1.1）----
    F.variant = map_variant_from_beta(F.beta);

    % ---- features ----
    gamma = 1/7;
    Fcore = compute_features_from_wave(W.f, W.init, W.beta, gamma, params);
    fn = fieldnames(Fcore);
    for i=1:numel(fn), F.(fn{i}) = Fcore.(fn{i}); end

catch ME
    wave_id = sprintf('%s-%02d', char(W.province), W.wave_idx);
    st = ME.stack;
    if ~isempty(st)
        warning('wave %s failed at %s:%d -> %s', wave_id, st(1).name, st(1).line, ME.message);
    else
        warning('wave %s failed: %s', wave_id, ME.message);
    end
    F = [];
end
end

% ---- helper: 用 beta 映射变异株（向下截到 1 位小数）----
function lab = map_variant_from_beta(beta)
b = double(beta);
q = floor((b + 1e-12) * 10) / 10;
if abs(q - 0.3) < 1e-6
    lab = "Alpha";
elseif abs(q - 0.7) < 1e-6
    lab = "Delta";
elseif abs(q - 1.1) < 1e-6
    lab = "Omicron";
else
    lab = "Unknown";
end
end
