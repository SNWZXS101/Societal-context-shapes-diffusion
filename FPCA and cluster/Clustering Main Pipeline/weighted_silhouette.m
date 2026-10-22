function s = weighted_silhouette(D, labels, w)
% D: [N x N], labels: [N x 1], w: weights
N = size(D,1);
w = w(:); w = w / sum(w);
s = nan(N,1);
for i=1:N
    li = labels(i);
    in_i = (labels==li);
    out  = (labels~=li);
    % a(i): avg intra-cluster distance
    if sum(in_i)>1
        ai = sum(D(i,in_i).*w(in_i)') / max(eps, sum(w(in_i))-w(i));
    else
        ai = NaN;
    end
    % b(i): min avg distance to other clusters
    bj = inf;
    for L = unique(labels(out))'
        idx = (labels==L);
        bj = min(bj, sum(D(i,idx).*w(idx)')/max(eps, sum(w(idx))));
    end
    s(i) = (bj - ai) / max(bj, ai + eps);
end
end
