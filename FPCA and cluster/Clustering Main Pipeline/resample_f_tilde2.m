
function [F2, ok] = resample_f_tilde2(W, u)
% f~2 = f / (β SI/N) ; need SI/N (or incidence/β)
% If waves_dataset contains incidence series 'inc' per wave, we can use inc = β SI/N => f~2 = f / inc
Nw = numel(W.f);
M = numel(u); F2 = nan(M,Nw); ok=false;
has_inc = isfield(W, 'inc'); % <- 如果有观测新增
if ~has_inc
    % 尝试从 W.Y 或 W.incidence 等字段兼容
    candidates = {'inc','incidence','new_cases','Y'};
    for c = candidates
        if isfield(W, c{1}), W.inc = W.(c{1}); has_inc = true; break; end
    end
end
if ~has_inc, return; end

for i=1:Nw
    f = double(W.f{i}(:));
    t = double(W.T{i}(:));
    inc = double(W.inc{i}(:));
    if numel(t)>=2 && range(t)>0 && numel(inc)==numel(t)
        t0 = (t - min(t)) / (max(t)-min(t));
        ft = f ./ max(inc, 1e-9); % f / (β SI/N)
        F2(:,i) = interp1(t0, ft, u, 'pchip','extrap');
        ok = true;
    end
end
end
