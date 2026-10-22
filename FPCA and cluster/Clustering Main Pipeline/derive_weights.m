function [omega, used_col] = derive_weights(Tfeat)
% 仅从 features 表读取权重：严格使用 w_tail_14
assert(istable(Tfeat), 'derive_weights_from_features expects a table.');
vn = string(Tfeat.Properties.VariableNames);
hit = find(strcmpi(vn, "w_tail_14"), 1);
assert(~isempty(hit), 'features.mat 中找不到权重列 w_tail_14');
omega = double(Tfeat.(hit));
omega = omega(:);
used_col = "w_tail_14";
% 诊断：非正或缺失
bad = ~isfinite(omega) | omega<0;
if any(bad)
    warning('Found %d invalid weights in w_tail_14; replaced with median.', sum(bad));
    med = median(omega(~bad));
    omega(bad) = med;
end
end
