function protos = dba_prototypes(F, labels, K, w, n_iter)
% F: [M x N] curves on same grid
protos = zeros(size(F,1), K);
for k=1:K
    idx = (labels==k);
    if ~any(idx)
        protos(:,k) = mean(F,2);
        continue;
    end
    curves = F(:, idx);
    wk = w(idx);
    protos(:,k) = dba_barycenter(curves, wk, n_iter);
end
end

