function [W, Tfeat] = load_from_features(features_mat)
% LOAD_FROM_FEATURES  从 features.mat 稳健读取 Features + f/init
% 输出结构：
%   W.f{i}      : 每波的 double 列向量 f_i(t)
%   W.T{i}      : 与 f_i 对齐的时间（这里用 1:len(f) 生成，主流程会归一化）
%   W.init{i}   : 每波的 init struct（含 Se0/S0/I0/R0 等）
%   W.province  : N×1 string，若缺失则填 "(Unknown)"
%   W.variant   : N×1 string，若缺失则填 "(Unknown)"
%   W.wave_id   : N×1 string，若缺失则用 province-wave_idx 生成
%   W.wave_idx  : N×1 double，若缺失则用 (1:N) 生成
%
% 依赖：features.mat 内包含变量：
%   - Features   (table)
%   - f_cell     (cell, 每单元一个 double 向量)
%   - init_cell  (cell, 每单元一个 struct)

S = load(features_mat);
if ~isfield(S,'Features')
    error('features.mat 缺少变量 "Features"（table）');
end
if ~isfield(S,'f_cell') || ~isfield(S,'init_cell')
    error(['features.mat 必须同时包含 "f_cell" 与 "init_cell"。\n' ...
           '请在特征导出阶段使用：save(outMat, ''Features'',''f_cell'',''init_cell'',''-v7.3'');']);
end

Tfeat = S.Features;
Fcell = S.f_cell;
Icell = S.init_cell;

% 基本校验
if ~iscell(Fcell) || ~iscell(Icell)
    error('"f_cell"/"init_cell" 必须为 cell。');
end
N = height(Tfeat);
if numel(Fcell) ~= N || numel(Icell) ~= N
    error('f_cell/init_cell 与 Features 行数不一致：%d vs %d', numel(Fcell), N);
end

% ---- 组装 W ----
W = struct();
W.f      = cell(N,1);
W.T      = cell(N,1);
W.init   = cell(N,1);

% 元信息列：宽松匹配列名（不强制存在）
prov = pick_string_col(Tfeat, ["province","Province","prov","region","Region"]);
varn = pick_string_col(Tfeat, ["variant","Variant","strain","Strain"]);
wid  = pick_string_col(Tfeat, ["wave_id","WaveID","id","ID","Wave_Id"]);
widx = pick_numeric_col(Tfeat, ["wave_idx","idx","index","Index","Wave_Idx"]);

if isempty(prov), prov = repmat(string("(Unknown)"), N, 1); end
if isempty(varn), varn = repmat(string("(Unknown)"), N, 1); end
if isempty(widx), widx = (1:N)'; end

% 若缺 wave_id，用 province-wave_idx 构造
if isempty(wid)
    wid = strings(N,1);
    for i=1:N
        wid(i) = sprintf('%s-%02d', to_char(prov(i)), round(widx(i)));
    end
end

% 填充每波的 f/T/init
for i = 1:N
    fi = Fcell{i};
    if ~(isnumeric(fi) && isvector(fi))
        error('第 %d 行的 f 不是数值向量。', i);
    end
    fi = double(fi(:));
    W.f{i} = fi;
    W.T{i} = (1:numel(fi))';  % 用索引当时间轴，后续会归一化到 [0,1]
    ini = Icell{i};
    if ~isstruct(ini)
        error('第 %d 行的 init 不是 struct。', i);
    end
    W.init{i} = ini;
end

% 元信息赋值（全部用 column vector）
W.province = string(prov(:));
W.variant  = string(varn(:));
W.wave_id  = string(wid(:));
W.wave_idx = double(widx(:));

end

% ---------- helpers ----------
function x = pick_string_col(T, candidates)
x = [];
vn = string(T.Properties.VariableNames);
for k = 1:numel(candidates)
    hit = find(strcmpi(vn, candidates(k)), 1);
    if ~isempty(hit)
        x = string(T.(hit));
        x = x(:);
        return;
    end
end
end

function x = pick_numeric_col(T, candidates)
x = [];
vn = string(T.Properties.VariableNames);
for k = 1:numel(candidates)
    hit = find(strcmpi(vn, candidates(k)), 1);
    if ~isempty(hit)
        x = double(T.(hit));
        x = x(:);
        return;
    end
end
end

function s = to_char(x)
if isstring(x)
    s = char(x);
else
    s = char(string(x));
end
end
