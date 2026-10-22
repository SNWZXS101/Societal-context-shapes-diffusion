function D = pairwise_dtw(F, varargin)
% PAIRWISE_DTW  成对 DTW 距离矩阵（带宽 + 可选并行，自动无并行回退）
% F: [M x N]
% Name-Value:
%   'Window'      : Sakoe–Chiba 带宽，默认 []
%   'UseParallel' : true/false，默认 true（若无并行工具箱则自动改成 false）

p = inputParser;
addParameter(p,'Window',[], @(z) isempty(z) || isnumeric(z));
addParameter(p,'UseParallel',true, @(z)islogical(z)||ismember(z,[0 1]));
parse(p, varargin{:});
win   = p.Results.Window;
dopar = logical(p.Results.UseParallel);

[~, N] = size(F);
D = zeros(N, N);

% —— 自动探测并行能力；没有就改为串行 ——
if dopar && ~hasPCT()
    warning('未检测到 Parallel Computing Toolbox，已改为串行计算。');
    dopar = false;
end

% 生成上三角索引
[I, J] = find(triu(true(N), 1));
P = numel(I);
dvals = zeros(P,1);

% 尝试开池（如已存在则不重复；失败则串行）
if dopar
    try
        pool = gcp('nocreate'); %#ok<NASGU>
        if isempty(gcp('nocreate')), parpool; end
    catch
        warning('无法创建并行池，改为串行计算。');
        dopar = false;
    end
end

if dopar
    parfor k = 1:P
        i = I(k); j = J(k);
        dvals(k) = dtw_distance(F(:,i), F(:,j), win);
    end
else
    for k = 1:P
        i = I(k); j = J(k);
        dvals(k) = dtw_distance(F(:,i), F(:,j), win);
    end
end

% 回填
for k = 1:P
    D(I(k), J(k)) = dvals(k);
    D(J(k), I(k)) = dvals(k);
end
D(1:N+1:end) = 0;
end

% ===== helpers =====
function tf = hasPCT()
% 尽量稳健地检测并行工具箱
tf = false;
try
    tf = ~isempty(ver('parallel')) && license('test','Distrib_Computing_Toolbox');
catch
    % 某些环境 license(...) 不可用，忽略
    tf = ~isempty(ver('parallel'));
end
end
