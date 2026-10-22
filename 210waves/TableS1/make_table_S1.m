function make_table_S1(varargin)
% 生成 Supplementary Table S1（Excel），汇总 210 波 PINN 结果
% 用法：
% make_table_S1('rootDir','.', 'featuresMat','features_aug.mat', 'outXlsx','Table_S1.xlsx');

%% --------- 参数 ---------
ip = inputParser;
addParameter(ip,'rootDir','.',@(s)ischar(s)||isstring(s));
addParameter(ip,'featuresMat','features_aug.mat',@(s)ischar(s)||isstring(s));
addParameter(ip,'outXlsx','Table_S1.xlsx',@(s)ischar(s)||isstring(s));
addParameter(ip,'smoothDays',5,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
parse(ip,varargin{:});
P = ip.Results;

rootDir  = char(P.rootDir);
outXlsx  = char(P.outXlsx);

%% --------- 读 Features（拿 wave_id / wave_no / variant / region / start_date 提示）---------
T = load(P.featuresMat);
if isfield(T,'Features'), F = T.Features;
else, error('features_mat 中未找到 Features 表'); end

need = {'province','wave_id','wave_no','variant','region','start_date'};

% 先把命中列名统一成标准名
for k = 1:numel(need)
    hit = find(strcmpi(F.Properties.VariableNames, need{k}), 1);
    if ~isempty(hit)
        F.Properties.VariableNames{hit} = need{k};
    end
end

% 再为缺失列按类型补一列
for k = 1:numel(need)
    if ~ismember(need{k}, F.Properties.VariableNames)
        switch need{k}
            case {'province','wave_id','variant','region'}
                F.(need{k}) = strings(height(F),1);
            case 'wave_no'
                F.(need{k}) = nan(height(F),1);
            case 'start_date'
                F.(need{k}) = NaT(height(F),1);
        end
    end
end

%% --------- 汇总表容器（不含 province 与 wave_no）---------
S = table();
S.wave_id          = strings(height(F),1);
S.variant_phase    = strings(height(F),1);
S.region           = strings(height(F),1);

S.start_date       = NaT(height(F),1);
S.end_date         = NaT(height(F),1);
S.duration_days    = nan(height(F),1);
S.obs_peak_date    = NaT(height(F),1);
S.obs_peak_cases   = nan(height(F),1);
S.obs_total_cases  = nan(height(F),1);

S.beta             = nan(height(F),1);
S.gamma            = nan(height(F),1);
S.init_Se0         = nan(height(F),1);
S.init_S0          = nan(height(F),1);
S.init_I0          = nan(height(F),1);
S.init_R0          = nan(height(F),1);
S.f_len            = nan(height(F),1);
S.f_mean           = nan(height(F),1);
S.f_max            = nan(height(F),1);

S.loss_total       = nan(height(F),1);
S.loss_data        = nan(height(F),1);
S.loss_physics     = nan(height(F),1);

% 来自 parameters.xlsx
S.iterations       = nan(height(F),1);
S.resume           = strings(height(F),1);
S.seed             = nan(height(F),1);

% 拟合误差
S.rmse_fit         = nan(height(F),1);
S.mae_fit          = nan(height(F),1);

%% --------- 主循环 ---------
for r = 1:height(F)
    prov   = string(F.province(r));
    wn     = double(F.wave_no(r));
    wavid  = string(F.wave_id(r));
    if ~isfinite(wn) || wn<=0, wn = parse_wave_no_from_id(wavid); end

    waveDir    = fullfile(rootDir, char(prov), sprintf('Results%02d', wn));
    resultsDir = fullfile(waveDir, 'results');

    % ---- daily_cases.xlsx ----
    dc = fullfile(waveDir, 'daily_cases.xlsx');
    date_obs = datetime.empty; y_obs = [];
    if isfile(dc)
        TT = readtable(dc,'VariableNamingRule','preserve');
        if size(TT,1)>=2
            date_obs = TT{2:end,1};
            if ~isdatetime(date_obs), date_obs = try_parse_datetime(date_obs); end
            y_obs = double(TT{2:end,2});
            if P.smoothDays>1
                y_obs = movmean(y_obs,[floor((P.smoothDays-1)/2) ceil((P.smoothDays-1)/2)],'omitnan');
            end
        end
    else
        warning('%s 缺少 daily_cases.xlsx', waveDir);
    end

    % ---- results 下文件 ----
    [tsFile, pmFile, lsFile] = locate_result_files(resultsDir);

    % ---- time_series：初值 + f(t) ----
    Se0=NaN; S0=NaN; I0=NaN; R0=NaN; f = [];
    if ~isempty(tsFile)
        TS = readtable(tsFile,'VariableNamingRule','preserve');
        vnames = string(TS.Properties.VariableNames);
        norm   = normalize_names(vnames);
        Se_col = find_first(norm, ["se","se0","se_0","expo","exposed"]);
        S_col  = find_first(norm, ["s","s0","s_0","sus","susceptible"]);
        I_col  = find_first(norm, ["i","i0","i_0","inf","infectious"]);
        R_col  = find_first(norm, ["r","r0","r_0","rec","removed","recovered"]);
        f_col  = find_first(norm, ["f","ft","f_t","phi","policy","u","control","f_t_"]); % 多给几种别名
        if ~isempty(Se_col), Se0 = safe_num(TS{1,Se_col}); end
        if ~isempty(S_col),  S0  = safe_num(TS{1,S_col});  end
        if ~isempty(I_col),  I0  = safe_num(TS{1,I_col});  end
        if ~isempty(R_col),  R0  = safe_num(TS{1,R_col});  end
        if ~isempty(f_col)
            f = double(TS{:,f_col});
        else
            % 兜底：挑一列像 f(t)（长、在 0~1.2）
            numCols = find(varfun(@isnumeric,TS,'OutputFormat','uniform'));
            for c = numCols
                vv = double(TS{:,c});
                if nnz(isfinite(vv))>=5
                    lo = prctile(vv,1,'all','omitnan'); hi = prctile(vv,99,'all','omitnan');
                    if lo>=-0.05 && hi<=1.20, f = vv; break; end
                end
            end
        end
    else
        warning('%s 缺少 time_series.xlsx（或匹配文件）', resultsDir);
    end

    % ---- parameters：beta/gamma + iterations/resume/seed ----
    beta = NaN; gamma = NaN; iterations = NaN; resumeStr = ""; seed = NaN;
    if ~isempty(pmFile)
        PP = readtable(pmFile,'VariableNamingRule','preserve');
        vn  = string(PP.Properties.VariableNames);
        nvn = normalize_names(vn);

        % 允许多种写法（含 β/γ）
        bcol = find_first(nvn, ["beta","b","β","par_beta","parbeta"]);
        gcol = find_first(nvn, ["gamma","g","γ","par_gamma","pargamma"]);
        if ~isempty(bcol) && height(PP)>=1, beta  = safe_num(PP{1,bcol}); end
        if ~isempty(gcol) && height(PP)>=1, gamma = safe_num(PP{1,gcol}); end

        itcol = find_first(nvn, ["iterations","epochs","maxepochs","max_iter","niter"]);
        rscol = find_first(nvn, ["resume","resumetraining","continue","cont"]);
        sdcol = find_first(nvn, ["seed","randomseed","rngseed"]);
        if ~isempty(itcol) && height(PP)>=1, iterations = safe_num(PP{1,itcol}); end
        if ~isempty(rscol) && height(PP)>=1, resumeStr  = as_bool_string(PP{1,rscol}); end
        if ~isempty(sdcol) && height(PP)>=1, seed       = safe_num(PP{1,sdcol}); end

        % Name/Value 两列式兜底
        if (isnan(iterations) || resumeStr=="" || ~isfinite(beta) || ~isfinite(gamma) || isnan(seed)) && width(PP)>=2
            names = normalize_text_col(PP{:,1}); vals = PP{:,2};
            if isnan(iterations)
                id = find(contains(names, ["iterations","epochs","maxepochs","max_iter","niter"]),1);
                if ~isempty(id), iterations = safe_num(vals(id)); end
            end
            if resumeStr==""
                id = find(contains(names, ["resume","resumetraining","continue","cont"]),1);
                if ~isempty(id), resumeStr = as_bool_string(vals(id)); end
            end
            if isnan(seed)
                id = find(contains(names, ["seed","randomseed","rngseed"]),1);
                if ~isempty(id), seed = safe_num(vals(id)); end
            end
            if ~isfinite(beta)
                id = find(contains(names, ["beta","par_beta","parbeta","β"]),1);
                if ~isempty(id), beta = safe_num(vals(id)); end
            end
            if ~isfinite(gamma)
                id = find(contains(names, ["gamma","par_gamma","pargamma","γ"]),1);
                if ~isempty(id), gamma = safe_num(vals(id)); end
            end
        end
    else
        warning('%s 缺少 parameters.xlsx（或匹配文件）', resultsDir);
    end

    % ---- losses：按最后一行取 total_loss / data_loss / phy_loss ----
    loss_total=NaN; loss_data=NaN; loss_physics=NaN;
    if ~isempty(lsFile)
        LL = readtable(lsFile,'VariableNamingRule','preserve');
        lvn = normalize_names(string(LL.Properties.VariableNames));
        % 先匹配你给的精确列名（规范化后：total_loss→totalloss，data_loss→dataloss，phy_loss→phyloss）
        tl = find_first(lvn, ["totalloss","total_loss","loss_total","total"]);
        dl = find_first(lvn, ["dataloss","data_loss","loss_data","data"]);
        pl = find_first(lvn, ["phyloss","phy_loss","physicsloss","loss_physics","physics","phys"]);
        if ~isempty(tl), loss_total   = safe_num(LL{end,tl}); end
        if ~isempty(dl), loss_data    = safe_num(LL{end,dl}); end
        if ~isempty(pl), loss_physics = safe_num(LL{end,pl}); end
    else
        warning('%s 缺少 losses.xlsx（或匹配文件）', resultsDir);
    end

    % ---- 观测统计 ----
    if ~isempty(date_obs)
        start_date = date_obs(1);
        end_date   = date_obs(end);
        duration   = days(end_date - start_date) + 1;
        [obs_peak_cases, idxpk] = max(y_obs);
        obs_peak_date = date_obs(min(idxpk,numel(date_obs)));
        obs_total     = nansum(y_obs);
    else
        start_date=NaT; end_date=NaT; duration=NaN; obs_peak_date=NaT; obs_peak_cases=NaN; obs_total=NaN;
    end

    % ---- f 统计 ----
    if ~isempty(f), f_mean = mean(f,'omitnan'); f_max = max(f,[],'omitnan'); f_len = numel(f);
    else, f_mean=NaN; f_max=NaN; f_len=NaN; end

    % ---- ODE 验证误差 ----
    rmse=NaN; mae=NaN;
    try
        if ~isempty(f) && all(isfinite([Se0 S0 I0 R0 beta gamma])) && ~isempty(date_obs)
            out = integrate_seir_odeint_style(f, struct('Se0',Se0,'S0',S0,'I0',I0,'R0',R0), beta, gamma);
            date0 = date_obs(1);
            date_ode = date0 + days(out.t(:));
            y_ode    = max(0, out.I_ode(:));
            [d_common, ia, ib] = intersect(date_obs, date_ode);
            if ~isempty(d_common)
                yA = y_obs(ia); yB = y_ode(ib);
                rmse = sqrt(mean((yA - yB).^2,'omitnan'));
                mae  = mean(abs(yA - yB),'omitnan');
            end
        end
    catch ME
        warning('ODE 校验失败：%s — %s', char(wavid), ME.message);
    end

    % ---- 写入 ----
    S.wave_id(r)         = wavid;
    S.variant_phase(r)   = string(F.variant(r));
    S.region(r)          = string(F.region(r));

    S.start_date(r)      = start_date;
    S.end_date(r)        = end_date;
    S.duration_days(r)   = duration;
    S.obs_peak_date(r)   = obs_peak_date;
    S.obs_peak_cases(r)  = obs_peak_cases;
    S.obs_total_cases(r) = obs_total;

    S.beta(r)            = beta;
    S.gamma(r)           = gamma;
    S.init_Se0(r)        = Se0;
    S.init_S0(r)         = S0;
    S.init_I0(r)         = I0;
    S.init_R0(r)         = R0;
    S.f_len(r)           = f_len;
    S.f_mean(r)          = f_mean;
    S.f_max(r)           = f_max;

    S.loss_total(r)      = loss_total;
    S.loss_data(r)       = loss_data;
    S.loss_physics(r)    = loss_physics;

    S.iterations(r)      = iterations;
    S.resume(r)          = resumeStr;
    S.seed(r)            = seed;

    S.rmse_fit(r)        = rmse;
    S.mae_fit(r)         = mae;
end

%% --------- 排序 & 日期格式 ---------
[provKey, wnoKey] = split_wave_id_for_sort(S.wave_id);
[~,ord] = sortrows([provKey, wnoKey]);
S = S(ord,:);

% 日期列统一去时分秒 + 格式
S.start_date    = dateshift(S.start_date,'start','day');
S.end_date      = dateshift(S.end_date,'start','day');
S.obs_peak_date = dateshift(S.obs_peak_date,'start','day');
S.start_date.Format    = 'yyyy-MM-dd';
S.end_date.Format      = 'yyyy-MM-dd';
S.obs_peak_date.Format = 'yyyy-MM-dd';

%% --------- 写 Excel（表注 + 数据） ---------
notes = {
'Table S1. Summary of PINN fits and ODE validation for 210 province–waves (2020–2022).';
'Columns include wave identity, variant phase and region labels; observed series statistics (start/end, duration, peak, total);';
'PINN/ODE quantities (beta, gamma, initial states, SMF f(t) length and summary), training losses (last record),';
'and validation errors (RMSE/MAE) against observed daily cases. Dates use yyyy-MM-dd.';
'Note: iterations and resume are read from parameters.xlsx. When resume = TRUE, iterations indicate the last epoch of the resumed run,';
'and do not necessarily reflect the cumulative training epochs across all runs. Seed is the random initialization seed.';
};

if isfile(outXlsx), delete(outXlsx); end
writecell(notes, outXlsx, 'Sheet','TableS1', 'Range','A1');                 % 1~6 行：表注
writetable(S,   outXlsx, 'Sheet','TableS1', 'Range','A6');                  % 6 行开始写表头+数据（不会覆盖表注）

fprintf('[Table S1] written to %s\n', outXlsx);
end

%% ================= helpers =================
function wn = parse_wave_no_from_id(wavid)
wn = NaN;
try
    s = char(wavid);
    k = regexp(s,'-(\d+)$','tokens','once');
    if ~isempty(k), wn = str2double(k{1}); end
catch
end
end

function [provKey, wnoKey] = split_wave_id_for_sort(wave_id)
wave_id = string(wave_id);
provKey = strings(size(wave_id));
wnoKey  = zeros(size(wave_id));
for i=1:numel(wave_id)
    tok = regexp(wave_id(i),'^([A-Za-z]+)-(\d+)$','tokens','once');
    if ~isempty(tok)
        provKey(i) = lower(string(tok{1}));
        wnoKey(i)  = str2double(tok{2});
    else
        provKey(i) = lower(wave_id(i)); wnoKey(i)=0;
    end
end
end

function dt = try_parse_datetime(x)
try, dt = datetime(x); return; catch, end
try, dt = datetime(x,'InputFormat','yyyy-M-d'); return; catch, end
try, dt = datetime(x,'InputFormat','yyyy-MM-dd'); return; catch, end
dt = datetime(x,'ConvertFrom','datenum');
end

function [tsFile, pmFile, lsFile] = locate_result_files(resultsDir)
tsFile=''; pmFile=''; lsFile='';
if ~isfolder(resultsDir), return; end
L = dir(fullfile(resultsDir,'*time*series*.xlsx'));
if isempty(L), L = dir(fullfile(resultsDir,'time_series.xlsx')); end
if ~isempty(L), [~,i]=max([L.datenum]); tsFile = fullfile(L(i).folder,L(i).name); end

L = dir(fullfile(resultsDir,'*param*.xlsx'));
if isempty(L), L = dir(fullfile(resultsDir,'parameters.xlsx')); end
if ~isempty(L), [~,i]=max([L.datenum]); pmFile = fullfile(L(i).folder,L(i).name); end

L = dir(fullfile(resultsDir,'*loss*.xlsx'));
if isempty(L), L = dir(fullfile(resultsDir,'losses.xlsx')); end
if ~isempty(L), [~,i]=max([L.datenum]); lsFile = fullfile(L(i).folder,L(i).name); end
end

function nn = normalize_names(vnames)
nn = lower(string(vnames));
nn = regexprep(nn,'[\s%()\[\]\{\}\-\.\\/]+','');  % 去空白/符号：total_loss -> totalloss
nn = strrep(nn,'β','beta'); nn = strrep(nn,'γ','gamma'); nn = strrep(nn,'φ','phi');
nn = strrep(nn,'susceptibles','susceptible');
nn = strrep(nn,'infected','infectious');
nn = strrep(nn,'removed','recovered');
end

function idx = find_first(names, keys)
idx = [];
keys = string(keys);
for k=1:numel(keys)
    j = find(names == lower(keys(k)), 1); % 先精确
    if ~isempty(j), idx=j; return; end
end
for k=1:numel(keys)
    j = find(contains(names, lower(keys(k))), 1); % 再包含
    if ~isempty(j), idx=j; return; end
end
end

function n = safe_num(x)
if iscell(x), x = x{1}; end
if isstring(x) || ischar(x), x = str2double(x); end
n = double(x);
end

function s = normalize_text_col(col)
if iscell(col) || isstring(col) || ischar(col)
    s = lower(string(col));
else
    s = lower(string(string(col)));
end
s = regexprep(s,'[\s%()\[\]\{\}\-\.\\/]+','');
s = strrep(s,'β','beta'); s = strrep(s,'γ','gamma'); s = strrep(s,'φ','phi');
end

function out = as_bool_string(val)
if islogical(val)
    out = string(val);
elseif isnumeric(val)
    out = string(val~=0);
elseif isstring(val) || ischar(val)
    s = strip(lower(string(val)));
    if any(strcmp(s,["true","t","yes","y","on","1"]))
        out = "true";
    elseif any(strcmp(s,["false","f","no","n","off","0"]))
        out = "false";
    else
        out = string(val);
    end
else
    out = string(val);
end
out = upper(out);
end
