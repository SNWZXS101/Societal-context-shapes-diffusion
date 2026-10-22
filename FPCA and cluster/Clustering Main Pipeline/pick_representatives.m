function repr = pick_representatives(F, labels, prototypes, Krep)
if nargin<4, Krep=6; end
K = size(prototypes,2);
repr = cell(K,1);
for k=1:K
    idx = find(labels==k);
    if isempty(idx), repr{k} = []; continue; end
    % distance to prototype
    d = arrayfun(@(j) norm(F(:, idx(j)) - prototypes(:,k)), 1:numel(idx));
    [~,ord] = sort(d,'ascend');
    take = idx(ord(1:min(Krep,numel(ord))));
    repr{k} = take(:)';
end
end