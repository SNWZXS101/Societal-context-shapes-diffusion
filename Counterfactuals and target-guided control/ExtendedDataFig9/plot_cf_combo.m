function plot_cf_combo(out_root, wave_id, varargin)
% PLOT_CF_COMBO
% 从 counterfactual_main 的输出直接读表/时序并画图（不做任何重算）
% - 上图：f(t)（基线 + 场景）
% - 下图：inc(t)=beta*S*I/N（基线 + 场景）
%
% 用法：
%   plot_cf_combo('results_cf', 'Anhui-01', ...
%       'Scenarios', {'TGC_STRONG','NEG_STRONG','TGC_LONG','ADV3'}, ...
%       'SavePNG', true);

% ---------- 参数 ----------
p = inputParser; p.FunctionName = mfilename;
addRequired(p,'out_root',@(x)ischar(x)||isstring(x));
addRequired(p,'wave_id', @(x)ischar(x)||isstring(x));
addParameter(p,'Scenarios',{},@(x)iscellstr(x)||all(string(x)==string(x)));
addParameter(p,'SavePNG',false,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'LineWidth',1.8,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'FontSize',11,@(x)isnumeric(x)&&isscalar(x)&&x>=6);
parse(p,out_root,wave_id,varargin{:});
A = p.Results;

out_root = char(out_root);
wave_id  = char(wave_id);
fig_dir  = fullfile(out_root,'figures');
if ~isfolder(fig_dir), mkdir(fig_dir); end

% ---------- 读取结果表 ----------
res_path = fullfile(out_root,'cf_results.mat');
if ~isfile(res_path)
    error('找不到结果文件：%s（请先运行 counterfactual_main）', res_path);
end
S = load(res_path,'Results','Scenarios');
Results = S.Results;

% ---------- 读取时序 ----------
ts_path = fullfile(out_root,'timeseries', sprintf('wave_%s.mat', sanitize(wave_id)));
if ~isfile(ts_path)
    error('找不到该波的时序：%s', ts_path);
end
TSS = load(ts_path, 'TS');
TS  = TSS.TS;

t   = TS.t(:);
f0  = TS.f(:);
inc0= TS.inc(:);

% ---------- 场景列表处理 ----------
scns = string(A.Scenarios);
if isempty(scns)
    % 若未指定，默认用结果表里该波出现过的场景（保持顺序）
    scns = unique(string(Results.scenario(strcmp(string(Results.wave_id), string(wave_id)))), 'stable');
end
scns = scns(:);

% ---------- 准备颜色 & 线型 ----------
clr = lines(max(1,numel(scns)));   % 不指定颜色也行，这里给个固定的 palette
lw0 = A.LineWidth;

% ---------- 找到该波的结果行（用于 legend 标注 %Avoided） ----------
mask_wave = strcmp(string(Results.wave_id), string(wave_id));
ResW = Results(mask_wave, :);

% ---------- 画图 ----------
tlo = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
set(gcf,'Color','w');

% ---- 上：f(t) ----
ax1 = nexttile(1); hold(ax1,'on');
plot(ax1, t, f0, 'k-', 'LineWidth', lw0+0.4, 'DisplayName','Baseline f');

leg_text = {}; leg_h = [];

if isfield(TS,'meta') && isfield(TS.meta,'t_cross') && ~isempty(TS.meta.t_cross) && ~isnan(TS.meta.t_cross)
    xline(ax1, TS.meta.t_cross, ':', 'Color',[0.4 0.4 0.4], ...
          'LineWidth',1.5, 'HandleVisibility','off');
end

if isfield(TS,'meta') && isfield(TS.meta,'t0') && ~isnan(TS.meta.t0)
    xline(ax1, TS.meta.t0, '--', 'Color',[0.85 0.2 0.2], 'LineWidth',1.5, 'HandleVisibility','off');
end
% t_cross 维持原样的灰色点划线


% 逐场景叠加 f
for i = 1:numel(scns)
    tag = upper(strtrim(scns(i)));
    f_field = "f_"+tag;
    if ~isfield(TS, f_field)
        warning('场景 %s：在 TS 中未找到字段 %s，已跳过。', tag, f_field);
        continue;
    end
    
    f_scn = TS.(f_field)(:);

    % 结果表中取 avoided_pct（若有）
    pct = lookup_pct(ResW, tag);
    if isnan(pct)
        lbl = sprintf('%s', tag);
    else
        lbl = sprintf('%s (%.1f%% avoided)', tag, pct);
    end

    h = plot(ax1, t, f_scn, '-', 'LineWidth', lw0, 'Color', clr(i,:), 'DisplayName', lbl);
    leg_h(end+1) = h; %#ok<AGROW>
end
xlabel(ax1,'Day'); ylabel(ax1,'f(t)');
title(ax1, sprintf('%s — SMF edits', wave_id));
grid(ax1,'on'); set(ax1,'FontSize',A.FontSize);
lg1 = legend(ax1,'Location','best'); lg1.Box = 'off';

% ---- 下：inc(t) ----

ax2 = nexttile(2); hold(ax2,'on');
plot(ax2, t, inc0, 'k-', 'LineWidth', lw0+0.4, 'DisplayName','Baseline inc');

for i = 1:numel(scns)
    tag = upper(strtrim(scns(i)));
    inc_field = "inc_"+tag;
    if ~isfield(TS, inc_field), continue; end
    inc_scn = TS.(inc_field)(:);
    % 下图也带名字，独立 legend
    plot(ax2, t, inc_scn, '-', 'LineWidth', lw0, 'Color', clr(i,:), 'DisplayName', char(tag));
end

xlabel(ax2,'Day'); ylabel(ax2,'inc(t) = \beta S I / N');
title(ax2, sprintf('%s — Incidence', wave_id));
grid(ax2,'on'); set(ax2,'FontSize',A.FontSize);
lg2 = legend(ax2,'Location','best'); lg2.Box = 'off';

% ---- 保存 ----
if A.SavePNG
    png_path = fullfile(fig_dir, sprintf('combo_%s.png', sanitize(wave_id)));
    exportgraphics(tlo, png_path, 'Resolution', 300);
    fprintf('[plot_cf_combo] Saved → %s\n', png_path);
end

end % ===== 主函数结束 =====

% ---------- 小工具 ----------
function s = sanitize(s), s = regexprep(char(s), '[^a-zA-Z0-9-_]', '_'); end

function pct = lookup_pct(T, tag)
% 从 Results 的该波子表里找到对应场景的 avoided_pct（若无则 NaN）
pct = NaN;
if isempty(T), return; end
mask = strcmpi(string(T.scenario), string(tag));
if any(mask)
    v = T.avoided_pct(mask);
    if ~isempty(v), pct = double(v(1)); end
end
end
