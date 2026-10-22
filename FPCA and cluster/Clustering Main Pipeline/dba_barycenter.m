
function c = dba_barycenter(curves, w, n_iter)
% curves: [M x n] on same grid; w: [n x 1]
% Simple DBA on discrete grid using DTW alignment
[M, n] = size(curves);
w = w(:); w = w / sum(w);
% init as weighted mean
c = curves * (w / sum(w));
for it=1:n_iter
    % accumulate aligned points
    num = zeros(M,1); den = zeros(M,1);
    for j=1:n
        cj = curves(:,j);
        [p,q] = dtw_path(c, cj); % alignment path (indices in c and cj)
        for t=1:numel(p)
            num(p(t)) = num(p(t)) + w(j) * cj(q(t));
            den(p(t)) = den(p(t)) + w(j);
        end
    end
    upd = num ./ max(den, 1e-9);
    upd(~isfinite(upd)) = c(~isfinite(upd));
    c = 0.7*c + 0.3*upd;  % damping for stability
end
end