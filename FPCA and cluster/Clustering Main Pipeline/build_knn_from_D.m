function knn_idx = build_knn_from_D(D, k)
N = size(D,1);
[~,ord] = sort(D,2,'ascend');
knn_idx = ord(:,2:(k+1));
end