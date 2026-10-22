
function v = vi(L1, L2)
% Variation of Information (Meila 2007)
L1 = L1(:); L2 = L2(:);
n = numel(L1);
u1 = unique(L1); u2 = unique(L2);
P = zeros(numel(u1), numel(u2));
for i=1:numel(u1)
    for j=1:numel(u2)
        P(i,j) = sum(L1==u1(i) & L2==u2(j)) / n;
    end
end
p1 = sum(P,2); p2 = sum(P,1);
H = @(p) -nansum(p .* log2(p + eps));
I = 0;
for i=1:numel(u1)
    for j=1:numel(u2)
        if P(i,j)>0
            I = I + P(i,j)*log2(P(i,j)/(p1(i)*p2(j)));
        end
    end
end
v = H(p1) + H(p2) - 2*I;
end