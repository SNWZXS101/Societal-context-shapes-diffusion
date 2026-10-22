
% smf_fpca_cluster_main('features_mat','features.mat', ...
%                       'M',200,'K_default',5,'K_grid',3:7, ...
%                       'DTW_Window',0.1,'UseParallel',false);

smf_fpca_cluster_main('features_mat','features.mat', ...
    'M',200,'K_grid',3:8,'K_default',3, ...
    'NormMode','Aplus', ...
    'DTW_Window',0.1,'UseParallel',true, ...
    'MultipartAsSeparate',true,'MultipartLabelAs','K+1');

% 
% S = load('smf_fpca_cluster_results.mat');
% [minA,maxA,meanA] = deal(min(S.Aplus_after), max(S.Aplus_after), mean(S.Aplus_after))
% 
% % 簇规模
% labels = S.cluster.labels; K = S.cluster.K;
% counts = arrayfun(@(k) sum(labels==k), 1:K)
% 
% % 如果还出现 {207,1,2} 这类高度不均衡，再看那几个样本的曲线是否“形状确实异类”：
% u = S.u; F = S.F_norm;
% abn = find(ismember(labels, find(counts<=2)));
% figure; hold on; 
% plot(u, F(:, abn), 'r-', 'LineWidth', 1.5); 
% plot(u, mean(F(:, setdiff(1:numel(labels),abn)), 2), 'k--', 'LineWidth', 1.2);
% legend('可疑单点/小簇曲线','其余均值'); grid on; title('异常 vs 均值（已在统一网格 A^+ 归一）');
