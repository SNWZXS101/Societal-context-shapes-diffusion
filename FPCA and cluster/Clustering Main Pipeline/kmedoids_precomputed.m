function [idx, medoid_idx] = kmedoids_precomputed(D, K)
% KMEDOIDS_PRECOMPUTED  K-medoids with precomputed distances (robust init).
N = size(D,1);
D = double(D);
D = 0.5*(D + D.');           % symmetrize
D(1:N+1:end) = 0;

% 尝试内置（若可用）
try
    opts = statset('MaxIter',200,'Display','off');
    [idx, medoid_idx] = kmedoids(D, K, 'Distance','precomputed', ...
                                 'Options',opts, 'Replicates',10);
    return;
catch
    % 回退到自实现 PAM
end

[idx, medoid_idx] = pam_precomputed(D, K, 10, 200);
end

% ================= PAM (precomputed) =================
function [idx, medoids] = pam_precomputed(D, K, replicates, maxiter)
if nargin<3, replicates=5; end
if nargin<4, maxiter=200; end
N = size(D,1);
best_cost = inf; best_idx = []; best_medoids = [];

for rep = 1:replicates
    % ---- 健壮的 k-medoids++ 初始化 ----
    med = zeros(K,1);
    chosen = false(N,1);
    med(1) = randi(N); chosen(med(1)) = true;
    minD = D(:, med(1));

    for k = 2:K
        p = minD.^2;
        p(chosen) = 0;                 % 已选样本不能再选
        s = sum(p);
        if ~(s > 0)
            % 退化情况：所有点距离已选集合都近似 0
            % 选“最远”的未选样本；若仍不分明，则随机选一个未选
            [~, cand] = max(minD + 1e-12*rand(N,1));
            if chosen(cand)
                avail = find(~chosen);
                cand = avail(randi(numel(avail)));
            end
        else
            p = p / s;
            cand = randsample(N, 1, true, p);
            % 浮点误差保护
            if chosen(cand)
                avail = find(~chosen);
                cand = avail(randi(numel(avail)));
            end
        end
        med(k) = cand;
        chosen(cand) = true;
        minD = min(minD, D(:, cand));
    end

    % ---- 迭代细化 ----
    [idx0, cost0] = assign_by_medoids(D, med);
    improved = true; it = 0;
    while improved && it < maxiter
        improved = false; it = it + 1;
        for k = 1:K
            members = find(idx0==k);
            if numel(members) <= 1, continue; end
            curm = med(k);
            bestm = curm; bestc = cluster_cost(D, members, curm);
            for cand = members(:)'
                if cand == curm, continue; end
                cst = cluster_cost(D, members, cand);
                if cst < bestc - 1e-12
                    bestc = cst; bestm = cand;
                end
            end
            if bestm ~= curm
                med(k) = bestm; improved = true;
            end
        end
        [idx0, cost0] = assign_by_medoids(D, med);
    end

    if cost0 < best_cost
        best_cost = cost0; best_idx = idx0; best_medoids = med;
    end
end

idx = best_idx; medoids = best_medoids;
end

function [idx, total_cost] = assign_by_medoids(D, med)
N = size(D,1); K = numel(med);
dist2med = zeros(N,K);
for k=1:K, dist2med(:,k) = D(:, med(k)); end
[~, idx] = min(dist2med, [], 2);
total_cost = 0;
for k=1:K
    members = (idx==k);
    if any(members), total_cost = total_cost + sum(dist2med(members,k)); end
end
end

function c = cluster_cost(D, members, med)
c = sum(D(members, med));
end
