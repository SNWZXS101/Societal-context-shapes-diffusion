
function g = gap_like(D, labels, B)
% crude "gap-like": compare within-cluster dispersion to label permutations
if nargin<3, B=50; end
N = size(D,1);
wc = within_dispersion(D, labels);
perm = zeros(B,1);
for b=1:B
    perm(b) = within_dispersion(D, labels(randperm(N)));
end
g = mean(log(perm+eps)) - log(wc+eps);
end

function w = within_dispersion(D, labels)
w = 0;
for L = unique(labels(:))'
    idx = find(labels==L);
    if numel(idx)>1
        w = w + sum(squareform(D(idx,idx),'tovector'));
    end
end
end
