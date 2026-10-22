function a = ari(L1, L2)
% Adjusted Rand Index
L1 = L1(:); L2 = L2(:);
u1 = unique(L1); u2 = unique(L2);
n = numel(L1);
N = zeros(numel(u1), numel(u2));
for i=1:numel(u1)
    for j=1:numel(u2)
        N(i,j) = sum(L1==u1(i) & L2==u2(j));
    end
end
ni = sum(N,2); nj = sum(N,1);
sumComb = @(v) sum(v.*(v-1)/2);
t1 = sumComb(N(:));
t2 = sumComb(ni) * sumComb(nj) / sumComb(n);
den = 0.5*(sumComb(ni)+sumComb(nj)) - t2;
if den==0, a = 1; else, a = (t1 - t2) / den; end
end


