function d = dtw_distance(x, y, window)
% DTW_DISTANCE  L2-DTW 距离，带可选 Sakoe–Chiba 带宽约束
% x, y : 列向量
% window : 带宽（可选）
%   - 若 0<window<1    → 按比例约束（如 0.1 = ±10% 主对角）
%   - 若 window>=1 整数 → 绝对带宽（索引偏移不超过 window）
%   - 若空/缺省        → 无约束（完全 DTW）
%
% 返回:
%   d : sqrt 累积代价（与之前版本一致，便于替换）

x = x(:); y = y(:);
nx = numel(x); ny = numel(y);

if nargin < 3 || isempty(window)
    w = max(nx, ny);   % effectively unconstrained
else
    if window < 1, w = ceil(window * max(nx, ny));
    else,           w = floor(window);
    end
end

% 代价矩阵按带宽计算（避免全矩阵内存）
C = inf(nx+1, ny+1); C(1,1) = 0;

for i = 1:nx
    jmin = max(1, i - w);
    jmax = min(ny, i + w);
    xi = x(i);
    for j = jmin:jmax
        cost = (xi - y(j))^2;
        C(i+1, j+1) = cost + min([C(i, j+1), C(i+1, j), C(i, j)]);
    end
end

d = sqrt(C(nx+1, ny+1));
end
