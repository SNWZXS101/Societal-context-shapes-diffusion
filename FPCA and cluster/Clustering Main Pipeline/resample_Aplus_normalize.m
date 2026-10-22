function [Fgrid, Aplus_after] = resample_Aplus_normalize(fcell, Tcell, u, wave_id)
% 先插值到统一网格 u，再在 u 上按 A+ 归一（∫max(f,0)du = 1）
% 若某波 A+≤0，直接 error 并列出 wave_id（或行号），避免带错
N = numel(fcell); M = numel(u);
Fgrid = nan(M, N);
Aplus_after = nan(N,1);
bad = false(N,1);

for i = 1:N
    f = double(fcell{i}(:));
    t = double(Tcell{i}(:));
    if numel(t)>=2 && (max(t)>min(t))
        t0 = (t - min(t)) / (max(t) - min(t));
        fi = interp1(t0, f, u, 'pchip','extrap');
    else
        fi = repmat(mean(f,'omitnan'), M, 1);
    end
    % 只用正相面积 A+ 来确定尺度，但不要把负相裁掉
    Ap = trapz(u, max(fi, 0));      % 用正半部分算 A+
    if ~(Ap > 0)
        bad(i) = true;
        Fgrid(:,i) = 0;
        Aplus_after(i) = 0;
        continue;
    end
    Fgrid(:,i)     = fi / Ap;   
    Aplus_after(i) = 1;

end

if any(bad)
    if nargin<4 || isempty(wave_id), wave_id = (1:N)'; end
    list = wave_id(bad);
    error('Aplus normalization failed (A+==0) for %d wave(s). wave_id: %s', ...
          numel(list), strjoin(string(list(:))', ', '));
end
end
