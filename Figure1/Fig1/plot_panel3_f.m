function plot_panel3_f(features_mat, wave_sel)
% PLOT_PANEL3_FCELL  从 features.mat 的 f_cell 选一条 f(t) 并画出来（Panel 3 底图）
% 用法：
%   plot_panel3_fcell;                        % 画第 1 条
%   plot_panel3_fcell('features.mat', 12);    % 画第 12 条
%   plot_panel3_fcell('features.mat', "list");% 列出前若干条标签
%   plot_panel3_fcell('features.mat', "Zhejiang-Delta-2021Q3"); % 按标签匹配
%
% 期望 MAT 内含：
%   - f_cell : N×1 cell，每个元素为一条数值向量 f(t)
%   - Features : 表（可选），用于组成标签（province/variant/ID 等）
%   - start_date_vec : datetime（可选）

if nargin < 1 || isempty(features_mat), features_mat = 'features.mat'; end
assert(exist(features_mat,'file')==2, '未找到文件：%s（pwd=%s）', features_mat, pwd);

S = load(features_mat);
assert(isfield(S,'f_cell') && iscell(S.f_cell), 'MAT 中未找到 f_cell (cell 数组)。');

f_cell = S.f_cell(:);
N = numel(f_cell);

% ---------- 构造标签 ----------
labels = strings(N,1);
if isfield(S,'Features') && istable(S.Features) && height(S.Features)==N
    T = S.Features;

    prov = pickcol(T, ["Province","Prov","province","prov","Region","region"]);
    var  = pickcol(T, ["Variant","variant","Strain","strain","Period","period"]);
    wid  = pickcol(T, ["WaveID","waveID","waveId","ID","id","wave_name","name"]);

    if isfield(S,'start_date_vec') && numel(S.start_date_vec)==N
        sd = string(S.start_date_vec);
    else
        sd = strings(N,1);
    end

    for i = 1:N
        parts = string.empty(0,1);
        if ~isempty(prov), parts(end+1) = strtrim(string(prov(i))); end
        if ~isempty(var),  parts(end+1) = strtrim(string(var(i)));  end
        if ~isempty(wid),  parts(end+1) = strtrim(string(wid(i)));  end
        if sd(i)~="",      parts(end+1) = sd(i);                    end
        if isempty(parts), labels(i) = "wave"+i;
        else               labels(i) = strjoin(parts(parts~=""), "-");
        end
    end
else
    if isfield(S,'start_date_vec') && numel(S.start_date_vec)==N
        labels = string(S.start_date_vec);
    else
        labels = "wave" + string(1:N);
    end
end

% ---------- 列表模式 ----------
if nargin >= 2 && isstring(wave_sel) && wave_sel=="list"
    K = min(N, 25);
    fprintf('可用标签（前 %d/%d 条）：\n', K, N);
    for i = 1:K, fprintf('  %3d: %s\n', i, labels(i)); end
    return
end

% ---------- 选择一条 ----------
if nargin < 2 || isempty(wave_sel)
    idx = 1;
elseif isnumeric(wave_sel)
    idx = wave_sel;
else
    key = string(wave_sel);
    idx = find(labels == key, 1, 'first');
    if isempty(idx), idx = find(contains(labels, key), 1, 'first'); end
    assert(~isempty(idx), '未找到标签匹配：%s', key);
end
assert(idx>=1 && idx<=N, '索引越界：1..%d（你给的是 %d）', N, idx);

% ---------- 取 f 并画 ----------
f = f_cell{idx};
assert(isnumeric(f) && isvector(f), 'f_cell{%d} 不是数值向量。', idx);
f = f(:);
t = linspace(0,1,numel(f)).';

figure('Color','w'); hold on; box on;
plot(t, zeros(size(t)), '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1); % 零线
hp = plot(t, f, 'Color', [24 119 242]/255, 'LineWidth', 2);           % f(t)
xlabel('Normalized time t'); ylabel('f(t)');
title(sprintf('SMF — %s', labels(idx)), 'Interpreter','none');
legend(hp, {'f(t)'}, 'Location','best'); hold off;

end

% ===== 辅助：大小写不敏感地获取表列，返回真实列名再取值 =====
function col = pickcol(T, candidates)
col = [];
vn = T.Properties.VariableNames;
for k = 1:numel(candidates)
    nm = char(candidates(k));
    hit = find(strcmpi(vn, nm), 1, 'first');  % 大小写不敏感匹配
    if ~isempty(hit)
        realname = vn{hit};                    % 用真实列名取值
        v = T.(realname);
        % 只返回“可拼标签”的常见类型
        if iscell(v) || isstring(v) || ischar(v) || isnumeric(v) || isdatetime(v) || islogical(v) || iscategorical(v)
            % 统一成 string 向量，但保持长度
            try
                col = string(v);
            catch
                % 某些类型（如 cell 的混合）转不了，用原值
                col = v;
            end
            return
        end
    end
end
end
