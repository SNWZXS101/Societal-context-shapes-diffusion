function meta = build_meta(W, omega, Ftab, varargin)
% BUILD_META  汇总元数据；同时保存原始权重/归一化权重
% 可选输入：
%   'omega_raw' : 原始从特征表读取的权重（未归一化）
%   'used_col'  : 字段名（如 'w_tail'）

p = inputParser;
addParameter(p,'omega_raw',[],@(x)isnumeric(x)||isempty(x));
addParameter(p,'used_col',"",@(x)isstring(x) || ischar(x));
parse(p, varargin{:});
omega_raw = p.Results.omega_raw;
used_col  = string(p.Results.used_col);

meta = struct();
meta.province = string_fallback(W.province);
meta.variant  = string_fallback(W.variant);
meta.wave_id  = W.wave_id(:);
meta.wave_idx = W.wave_idx(:);

% 权重：原始 & 归一化
meta.weight        = omega(:);                  % 未归一化（直接供你核对）
if ~isempty(omega_raw)
    meta.weight_raw   = omega_raw(:);
else
    meta.weight_raw   = meta.weight;
end
meta.weight_source = used_col;                  % 记录使用了哪个列名（如 'w_tail'）

% 省内波序（若有）
meta.wave_order = [];
if istable(Ftab)
    vn = lower(string(Ftab.Properties.VariableNames));
    if any(vn=="wave_order")
        meta.wave_order = Ftab{:, find(vn=="wave_order",1)};
    elseif any(vn=="order")
        meta.wave_order = Ftab{:, find(vn=="order",1)};
    elseif any(vn=="wave_seq")
        meta.wave_order = Ftab{:, find(vn=="wave_seq",1)};
    end
end
end

function s = string_fallback(c)
if iscell(c), s = string(c(:));
elseif isstring(c), s = c(:);
else, s = strings(0,1);
end
end
