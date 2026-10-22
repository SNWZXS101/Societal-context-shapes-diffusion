function [Phi, scores, expl_var, mu] = weighted_fpca(F, w, ncomp)
% WEIGHTED_FPCA  Weighted functional PCA on curves matrix F [M x N]
% Inputs:
%   F     : [M x N] (columns = waves)
%   w     : [N x 1] sample weights (will be normalized inside)
%   ncomp : #PCs (default: min(10,N))
%
% Outputs:
%   Phi      : [M x r] eigenfunctions (PC loadings on the grid)
%   scores   : [N x r] PC scores per wave
%   expl_var : [r x 1] explained variance ratio
%   mu       : [M x 1] weighted mean curve

    % ---------- sanity & coercion ----------
    [M,N] = size(F);
    if nargin < 3 || isempty(ncomp), ncomp = min(10, N); end

    % Force w to be numeric column of length N
    if istable(w), w = table2array(w); end
    if iscell(w),  w = cellfun(@double, w); end
    w = double(w(:));
    if numel(w) ~= N
        error('weighted_fpca:weight_size', ...
              'Weight length mismatch: numel(w)=%d, N=%d', numel(w), N);
    end
    w(~isfinite(w) | w<=0) = 1e-6;
    w = w / sum(w);                 % normalize weights to sum=1

    % Replace NaNs in F by column means
    colmean = mean(F,'omitnan');
    for j = 1:N
        bad = ~isfinite(F(:,j));
        if any(bad), F(bad,j) = colmean(j); end
    end

    % ---------- weighted mean & centering ----------
    mu = F * w;                     % [M x 1]
    X  = F - mu;                    % [M x N], center columns

    % ---------- weighted covariance operator ----------
    % C = X * diag(w) * X'
    Xw = X .* (ones(M,1) * w');     % scale columns of X by weights
    C  = X * Xw';                   % [M x M]
    Cs = (C + C')/2;                % symmetrize for numeric stability

    % ---------- eigendecomposition ----------
    r  = min([ncomp, M, N]);
    try
        % Try eigs on symmetric matrix
        [U,S] = eigs(Cs, r, 'largestreal');
        lam = real(diag(S));
    catch
        % Fall back to full eig, then sort
        [U_full, S_full] = eig(Cs, 'vector');   % 'vector' returns eigenvalues as vector
        [lam_sorted, ord] = sort(real(S_full), 'descend');
        U = U_full(:, ord(1:r));
        lam = lam_sorted(1:r);
    end
    lam(lam<0) = 0;                 % numeric guard
    total = sum(lam) + eps;
    expl_var = lam / total;
    Phi = U;                        % [M x r]

    % ---------- scores ----------
    % L2 inner product on uniform grid (spacing absorbed in scaling)
    scores = (Phi' * X).';          % [N x r]
end
