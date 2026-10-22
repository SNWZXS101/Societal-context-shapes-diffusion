function [F_fm, fmax_after] = resample_fmax_normalize(Fcell, Tcell, u, wave_id)
% RESAMPLE_FMAX_NORMALIZE
% - 将每条 f(t) 插值到公共网格 u∈[0,1]
% - 用每条曲线的 f_max 做幅度归一（保留正负号）
%
% 输入:
%   Fcell : {N×1}，每单元一个 double 向量 f_i
%   Tcell : {N×1}，每单元对应时间向量 Ti（若为空或长度<2，则用等距索引）
%   u     : M×1 公共网格（0..1）
%   wave_id : (可选) {N×1} string/cellstr，用于告警信息
%
% 输出:
%   F_fm       : M×N 归一化后的矩阵
%   fmax_after : N×1，本条曲线用于归一的 f_max（用于日志/诊断）

N = numel(Fcell); M = numel(u);
F_fm = nan(M, N);
fmax_after = zeros(N,1);

if nargin < 4, wave_id = repmat({''}, N, 1); end
if isstring(wave_id), wave_id = cellstr(wave_id); end

for i = 1:N
    f = double(Fcell{i}(:));
    if isempty(f) || ~any(isfinite(f))
        F_fm(:,i) = 0; fmax_after(i) = 0; continue;
    end
    % 时间轴
    if ~isempty(Tcell) && numel(Tcell)>=i && ~isempty(Tcell{i})
        t = double(Tcell{i}(:));
        if numel(t) >= 2 && max(t) > min(t)
            t0 = (t - min(t)) / (max(t) - min(t));
        else
            t0 = linspace(0,1,numel(f))';
        end
    else
        t0 = linspace(0,1,numel(f))';
    end
    % 插值到公共网格
    fi = interp1(t0, f, u, 'pchip','extrap');

    % f_max（优先用正峰；若不存在正峰则用 |f| 的最大值兜底）
    fmax_pos = max(fi);
    if fmax_pos > 0 && isfinite(fmax_pos)
        s = fmax_pos;
    else
        s = max(abs(fi));
    end

    if ~(s > 0 && isfinite(s))
        % 全零或极端异常，直接清零
        F_fm(:,i) = 0; fmax_after(i) = 0;
        if ~isempty(wave_id{i})
            warning('FmaxNorm:ZeroScale', 'wave %s has non-positive scale; set to 0.', wave_id{i});
        end
    else
        F_fm(:,i) = fi / s;      % ✨ 保留负相，按 f_max 归一
        fmax_after(i) = s;
    end
end
end
