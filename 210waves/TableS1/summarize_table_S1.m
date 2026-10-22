function summarize_table_S1(xlsxPath, varargin)
% summarize_table_S1('Table_S1.xlsx', 'ExcludeMega', true, 'Decimals', 1)
% 读取 Table_S1.xlsx，输出 RMSE/MAE 的总体与分期统计（英文句子，可直接粘贴正文）。
%
% 参数：
%   xlsxPath      : Table_S1.xlsx 路径
% Name-Value：
%   'ExcludeMega' : true/false，是否剔除 Hubei-01 / Shanghai-02 / Jilin-03（默认 true）
%   'Decimals'    : 小数位数（默认 1）
%   'Quiet'       : 不打印分布个数补充（默认 false）

ip = inputParser;
addParameter(ip,'ExcludeMega',true,@islogical);
addParameter(ip,'Decimals',1,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(ip,'Quiet',false,@islogical);
parse(ip,varargin{:});
P = ip.Results;

% ----- 读取 -----
T = readtable(xlsxPath, 'Sheet','TableS1');  % A6 起是数据，readtable会自动识别
need = {'wave_id','variant_phase','rmse_fit','mae_fit'};
for k = 1:numel(need)
    if ~ismember(need{k}, T.Properties.VariableNames)
        error('Table_S1.xlsx 缺少列：%s', need{k});
    end
end

% 类型规范
T.variant_phase = string(T.variant_phase);
T.wave_id       = string(T.wave_id);
RMSE = double(T.rmse_fit);
MAE  = double(T.mae_fit);

% 剔除 NaN
valid = isfinite(RMSE) & isfinite(MAE);
T = T(valid,:); RMSE = RMSE(valid); MAE = MAE(valid);

% 是否剔除三大波
if P.ExcludeMega
    isMega = ismember(lower(T.wave_id), lower(["Hubei-01","Shanghai-02","Jilin-03"]));
    T = T(~isMega,:); RMSE = RMSE(~isMega); MAE = MAE(~isMega);
end

% 归一化分期名称
V = lower(strrep(T.variant_phase," ",""));
V(contains(V,["original","alpha"])) = "Original/Alpha";
V(contains(V,"delta"))              = "Delta";
V(contains(V,"omicron"))            = "Omicron";

% ----- 统计函数 -----
dq = @(x,p) quantile(x,p);
fmt = @(x) sprintf(['%0.',num2str(P.Decimals),'f'], x);
fmt_iqr = @(x) sprintf(['%s [%s–%s]'], fmt(median(x,'omitnan')), fmt(dq(x,0.25)), fmt(dq(x,0.75)));
fmt_mean_sd = @(x) sprintf(['%s ± %s'], fmt(mean(x,'omitnan')), fmt(std(x,'omitnan')));

% ----- 总体 -----
sent_overall = sprintf( ...
    ['Across 210 waves, reconstruction errors were modest (overall; ' ...
     'RMSE median [IQR] = %s cases/day, MAE median [IQR] = %s cases/day; ' ...
     'mean ± SD: RMSE %s, MAE %s; Table S1).'], ...
     fmt_iqr(RMSE), fmt_iqr(MAE), fmt_mean_sd(RMSE), fmt_mean_sd(MAE));

% ----- 分期（Original/Alpha, Delta, Omicron）-----
labels = ["Original/Alpha","Delta","Omicron"];
parts = strings(1, numel(labels));
for i = 1:numel(labels)
    mask = (V == labels(i));
    if any(mask)
        parts(i) = sprintf('%s: RMSE %s; MAE %s', ...
            labels(i), fmt_iqr(RMSE(mask)), fmt_iqr(MAE(mask)));
    else
        parts(i) = sprintf('%s: n=0', labels(i));
    end
end
sent_byperiod = sprintf('By variant period, %s; see Table S1.', strjoin(parts, '; '));

% ----- 可选数量补充 -----
if ~P.Quiet
    n_overall = height(T);
    n_oa = sum(V=="Original/Alpha"); n_de = sum(V=="Delta"); n_om = sum(V=="Omicron");
    fprintf('[n] overall=%d (Original/Alpha=%d, Delta=%d, Omicron=%d)\n', n_overall, n_oa, n_de, n_om);
    if P.ExcludeMega, fprintf('[note] three mega-waves excluded.\n'); end
end

% ----- 打印可直接粘贴的句子 -----
fprintf('\n%s\n', sent_overall);
fprintf('%s\n\n', sent_byperiod);

% 同时返回字符串（如果作为函数调用）
if nargout>0
    varargout{1} = sent_overall;
    if nargout>1, varargout{2} = sent_byperiod; end
end
end
