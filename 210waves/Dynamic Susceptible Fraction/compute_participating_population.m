function Summary = compute_participating_population(rootDir)
% COMPUTE_PARTICIPATING_POPULATION
%
% 统计每个疫情波次中：
%
%   1. 总人口池
%          P = Se + S + I + R
%
%   2. 初始传播核心人口
%          Core0 = S0 + I0 + R0
%
%   3. 最大累计动员人口
%          MobilizedMax = Se0 - min(Se)
%
%   4. 最大传播核心规模
%          CoreMax = Core0 + MobilizedMax
%                  = max(S + I + R)
%
%   5. 三个主要比例：
%          InitialCorePct
%          MobilizedPct
%          MaxCorePct
%
% 同时进行质量控制：
%   - Se+S+I+R 的人口守恒误差
%   - 两种 CoreMax 算法的一致性
%
% 使用：
%   Summary = compute_participating_population(pwd);
%
% 若本文件就在30省文件夹同级：
%   Summary = compute_participating_population();
%
% 输出：
%   participating_population_by_wave.xlsx
%   participating_population_overall.txt
%   participating_population_summary.mat

if nargin < 1 || isempty(rootDir)
    rootDir = pwd;
end

rootDir = char(rootDir);

%% --------------------------------------------------------
% 查找省份文件夹
% ---------------------------------------------------------

dd = dir(rootDir);
dd = dd([dd.isdir]);

excludeNames = {'.','..','figures','results','output','outputs'};
keep = ~ismember({dd.name}, excludeNames);
dd = dd(keep);

rows = {};

fprintf('\n===============================================\n');
fprintf('Scanning province folders...\n');
fprintf('Root directory: %s\n', rootDir);
fprintf('===============================================\n\n');

%% --------------------------------------------------------
% 遍历省份
% ---------------------------------------------------------

for ip = 1:numel(dd)

    province = dd(ip).name;
    provinceDir = fullfile(rootDir, province);

    % 只寻找 ResultsXX
    Rdirs = dir(fullfile(provinceDir, 'Results*'));
    Rdirs = Rdirs([Rdirs.isdir]);

    if isempty(Rdirs)
        continue;
    end

    % Results01, Results02 ... 按数字排序
    resultNums = nan(numel(Rdirs),1);

    for j = 1:numel(Rdirs)
        tok = regexp(Rdirs(j).name, ...
            '^Results(\d+)$', 'tokens', 'once');

        if ~isempty(tok)
            resultNums(j) = str2double(tok{1});
        end
    end

    valid = isfinite(resultNums);

    Rdirs = Rdirs(valid);
    resultNums = resultNums(valid);

    [resultNums, ord] = sort(resultNums);
    Rdirs = Rdirs(ord);

    fprintf('%-20s : %d wave folders\n', ...
        province, numel(Rdirs));

    %% ----------------------------------------------------
    % 每一波
    % -----------------------------------------------------

    for iw = 1:numel(Rdirs)

        resultName = Rdirs(iw).name;
        waveIdx = resultNums(iw);

        resultDir = fullfile(provinceDir, resultName);

        % 优先：
        % ResultsXX/results/time_series.xlsx
        tsFile = fullfile(resultDir, ...
            'results', 'time_series.xlsx');

        % 兼容：
        % ResultsXX/time_series.xlsx
        if ~isfile(tsFile)
            tsFile = fullfile(resultDir, ...
                'time_series.xlsx');
        end

        if ~isfile(tsFile)
            warning('SMF:MissingTimeSeries', ...
                'Missing time_series.xlsx: %s', resultDir);
            continue;
        end

        %% ------------------------------------------------
        % 读时序
        % -------------------------------------------------

        T = readtable(tsFile, ...
            'VariableNamingRule','preserve');

        names = string(T.Properties.VariableNames);

        cSe = find_column(names, "Se");
        cS  = find_column(names, "S");
        cI  = find_column(names, "I");
        cR  = find_column(names, "R");

        if any([cSe,cS,cI,cR] == 0)

            warning('SMF:MissingCompartments', ...
                ['Missing Se/S/I/R columns in %s.\n' ...
                 'Columns found: %s'], ...
                tsFile, strjoin(names, ', '));

            continue;
        end

        Se = double(T{:,cSe});
        S  = double(T{:,cS});
        I  = double(T{:,cI});
        R  = double(T{:,cR});

        Se = Se(:);
        S  = S(:);
        I  = I(:);
        R  = R(:);

        nT = numel(Se);

        if any([numel(S),numel(I),numel(R)] ~= nT)
            warning('SMF:LengthMismatch', ...
                'Length mismatch: %s', tsFile);
            continue;
        end

        %% ------------------------------------------------
        % 清理 NaN
        % -------------------------------------------------

        validRow = ...
            isfinite(Se) & ...
            isfinite(S)  & ...
            isfinite(I)  & ...
            isfinite(R);

        Se = Se(validRow);
        S  = S(validRow);
        I  = I(validRow);
        R  = R(validRow);

        if isempty(Se)
            continue;
        end

        %% =================================================
        % 核心统计
        % ==================================================

        % ---- 初值 ----

        Se0 = Se(1);
        S0  = S(1);
        I0  = I(1);
        R0  = R(1);

        % 总人口池
        P0 = Se0 + S0 + I0 + R0;

        if ~(isfinite(P0) && P0 > 0)
            warning('SMF:InvalidPopulation', ...
                'Invalid population for %s/%s', ...
                province, resultName);
            continue;
        end

        % --------------------------------------------
        % 总人口守恒检查
        % --------------------------------------------

        P_t = Se + S + I + R;

        conservationAbs = max(abs(P_t - P0));

        conservationRel = ...
            conservationAbs / P0;

        % --------------------------------------------
        % 初始传播核心
        % --------------------------------------------

        Core0 = S0 + I0 + R0;

        InitialCoreFrac = Core0 / P0;

        % --------------------------------------------
        % Se最低点
        % --------------------------------------------

        [SeMin, idxSeMin] = min(Se);

        % 累计净动员的最大值
        MobilizedMax = Se0 - SeMin;

        % 数值噪声保护
        if MobilizedMax < 0 && ...
                abs(MobilizedMax/P0) < 1e-10

            MobilizedMax = 0;
        end

        MobilizedFrac = MobilizedMax / P0;

        % --------------------------------------------
        % 方法 A：
        % 根据 Se 减少量得到 CoreMax
        % --------------------------------------------

        CoreMax_from_Se = ...
            Core0 + MobilizedMax;

        % --------------------------------------------
        % 方法 B：
        % 直接计算 S+I+R
        % --------------------------------------------

        Core_t = S + I + R;

        [CoreMax_direct, idxCoreMax] = ...
            max(Core_t);

        MaxCoreFrac = ...
            CoreMax_direct / P0;

        % 两种方法的一致性
        CoreConsistencyAbs = ...
            abs(CoreMax_direct - ...
                CoreMax_from_Se);

        CoreConsistencyRel = ...
            CoreConsistencyAbs / P0;

        % --------------------------------------------
        % 另外一个辅助指标：
        %
        % 当前仍可直接进入感染过程的人：
        % S + I
        %
        % R 已退出感染过程，但曾经属于传播核心
        % --------------------------------------------

        TransmissionActive_t = S + I;

        ActiveMax = max(TransmissionActive_t);

        ActiveMaxFrac = ...
            ActiveMax / P0;

        %% ------------------------------------------------
        % 时间位置
        % -------------------------------------------------

        if any(strcmpi(names,'Day'))
            day = double(T{validRow, ...
                find(strcmpi(names,'Day'),1)});

            day = day(:);

            tSeMin   = day(idxSeMin);
            tCoreMax = day(idxCoreMax);

        else
            tSeMin   = idxSeMin - 1;
            tCoreMax = idxCoreMax - 1;
        end

        %% ------------------------------------------------
        % 保存一行
        % -------------------------------------------------

        rows(end+1,:) = { ...
            string(province), ...
            waveIdx, ...
            string(resultName), ...
            P0, ...
            Se0, ...
            S0, ...
            I0, ...
            R0, ...
            Core0, ...
            InitialCoreFrac * 100, ...
            SeMin, ...
            MobilizedMax, ...
            MobilizedFrac * 100, ...
            CoreMax_from_Se, ...
            CoreMax_direct, ...
            MaxCoreFrac * 100, ...
            ActiveMax, ...
            ActiveMaxFrac * 100, ...
            tSeMin, ...
            tCoreMax, ...
            conservationRel, ...
            CoreConsistencyRel, ...
            string(tsFile) ...
            }; %#ok<AGROW>

    end
end

%% ========================================================
% 转为 table
% =========================================================

varNames = { ...
    'province', ...
    'wave_idx', ...
    'result_folder', ...
    'P_total', ...
    'Se0', ...
    'S0', ...
    'I0', ...
    'R0', ...
    'Core0', ...
    'InitialCorePct', ...
    'Se_min', ...
    'MobilizedMax', ...
    'MobilizedPct', ...
    'CoreMax_from_Se', ...
    'CoreMax_direct', ...
    'MaxCorePct', ...
    'ActiveSI_Max', ...
    'ActiveSI_MaxPct', ...
    't_Se_min', ...
    't_Core_max', ...
    'population_conservation_rel_error', ...
    'core_identity_rel_error', ...
    'source_file' ...
    };

Summary = cell2table(rows, ...
    'VariableNames', varNames);

%% ========================================================
% 汇总
% =========================================================

Nwave = height(Summary);

fprintf('\n===============================================\n');
fprintf('TOTAL WAVES FOUND = %d\n', Nwave);
fprintf('===============================================\n');

if Nwave == 0
    error('No valid wave data found.');
end

x = Summary.MaxCorePct;

fprintf('\nMaximum transmission-core population fraction:\n');
fprintf('  MIN       = %.4f %%\n', min(x));
fprintf('  P2.5      = %.4f %%\n', prctile(x,2.5));
fprintf('  P25       = %.4f %%\n', prctile(x,25));
fprintf('  MEDIAN    = %.4f %%\n', median(x));
fprintf('  P75       = %.4f %%\n', prctile(x,75));
fprintf('  P97.5     = %.4f %%\n', prctile(x,97.5));
fprintf('  MAX       = %.4f %%\n', max(x));

fprintf('\nInitial core population fraction:\n');
fprintf('  range     = %.4f %% -- %.4f %%\n', ...
    min(Summary.InitialCorePct), ...
    max(Summary.InitialCorePct));

fprintf('\nAdditional population mobilized during wave:\n');
fprintf('  range     = %.4f %% -- %.4f %%\n', ...
    min(Summary.MobilizedPct), ...
    max(Summary.MobilizedPct));

%% ========================================================
% 找极端波
% =========================================================

[~,imin] = min(Summary.MaxCorePct);
[~,imax] = max(Summary.MaxCorePct);

fprintf('\nWave with LOWEST maximum core fraction:\n');
disp(Summary(imin, ...
    {'province','wave_idx', ...
     'MaxCorePct','InitialCorePct','MobilizedPct'}));

fprintf('Wave with HIGHEST maximum core fraction:\n');
disp(Summary(imax, ...
    {'province','wave_idx', ...
     'MaxCorePct','InitialCorePct','MobilizedPct'}));

%% ========================================================
% QC
% =========================================================

fprintf('\nQuality control:\n');

fprintf('  Max population-conservation relative error = %.3e\n', ...
    max(Summary.population_conservation_rel_error));

fprintf('  Max core-identity relative error           = %.3e\n', ...
    max(Summary.core_identity_rel_error));

%% ========================================================
% 保存
% =========================================================

xlsxFile = fullfile(rootDir, ...
    'participating_population_by_wave.xlsx');

writetable(Summary, xlsxFile);

save(fullfile(rootDir, ...
    'participating_population_summary.mat'), ...
    'Summary');

txtFile = fullfile(rootDir, ...
    'participating_population_overall.txt');

fid = fopen(txtFile,'w');

fprintf(fid, ...
    'Number of valid waves: %d\n\n', Nwave);

fprintf(fid, ...
    'Maximum transmission-core population fraction:\n');

fprintf(fid, ...
    'Min: %.6f %%\n', min(x));

fprintf(fid, ...
    'Median: %.6f %%\n', median(x));

fprintf(fid, ...
    'Max: %.6f %%\n\n', max(x));

fprintf(fid, ...
    'Initial core fraction range: %.6f -- %.6f %%\n', ...
    min(Summary.InitialCorePct), ...
    max(Summary.InitialCorePct));

fprintf(fid, ...
    'Additional mobilized fraction range: %.6f -- %.6f %%\n', ...
    min(Summary.MobilizedPct), ...
    max(Summary.MobilizedPct));

fclose(fid);

fprintf('\nSaved:\n');
fprintf('  %s\n', xlsxFile);
fprintf('  %s\n', txtFile);
fprintf('  participating_population_summary.mat\n');

end

%% =========================================================
% helper
% ==========================================================

function idx = find_column(names, target)

idx = find(strcmpi(strtrim(names), target), ...
    1, 'first');

if isempty(idx)
    idx = 0;
end

end