function [F1, ok] = resample_f_tilde1(W, u)
% f~1 = (beta/N)*f; try to get N from init if present
Nw = numel(W.f);
M = numel(u); F1 = nan(M,Nw); ok=false;
for i=1:Nw
    f = double(W.f{i}(:));
    t = double(W.T{i}(:));
    beta = W.beta(i);
    Ntot = try_get_N(W, i);  % from init or []
    if isempty(Ntot) || ~isfinite(Ntot) || Ntot<=0
        % fall back to 1 (scales cancel in normalization across waves anyway)
        Ntot = 1;
    end
    ft = (beta / Ntot) * f;
    if numel(t)>=2 && range(t)>0
        t0 = (t - min(t)) / (max(t)-min(t));
        F1(:,i) = interp1(t0, ft, u, 'pchip','extrap');
        ok = true;
    end
end
end

