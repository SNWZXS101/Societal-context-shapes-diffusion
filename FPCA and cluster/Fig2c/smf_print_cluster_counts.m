function counts = smf_print_cluster_counts(labels, varargin)
% 打印聚类计数，便于快速核对是否一致
% 用法：
%   counts = smf_print_cluster_counts(labels)                  % 基本
%   counts = smf_print_cluster_counts(labels, 'Tag','Fig2a')   % 打标签
%   counts = smf_print_cluster_counts(labels, 'ExpectK',4)     % 期望 K（可选）
%
% 输出：
%   counts: 1xK 各簇计数（K = max(labels)）

p = inputParser;
p.addParameter('Tag','',@(x)ischar(x)||isstring(x));
p.addParameter('ExpectK',[],@(x)isempty(x)||(isscalar(x)&&x==round(x)));
p.parse(varargin{:});
tag = char(p.Results.Tag);
Kexp = p.Results.ExpectK;

labels = labels(:);
if any(~isfinite(labels))
    warning('labels 中存在 NaN/Inf，将被忽略计数。');
end

% 取实际 K（不要信 S.cluster.K）
L = labels(isfinite(labels));
if isempty(L)
    fprintf('[%s] labels 为空。\n', tag); counts = []; return;
end
K = max(L);
% 计数
counts = accumarray(max(1,min(K, round(L))), 1, [K 1])';  % 容错到 1..K
N = numel(L);

% 打印
if ~isempty(tag), tinfo = [' [' tag ']']; else, tinfo = ''; end
fprintf('%s Cluster counts (K=%d, N=%d): %s\n', tinfo, K, N, mat2str(counts));

% 细化校验
if ~isempty(Kexp) && K ~= Kexp
    warning('%s 实际 K=%d 与期望 K=%d 不一致。', tinfo, K, Kexp);
end
if sum(counts) ~= N
    warning('%s 计数之和(%d) 与 N(%d) 不一致。', tinfo, sum(counts), N);
end
% 每类逐条显示
for k = 1:K
    fprintf('  C%d: %d\n', k, counts(k));
end
end
