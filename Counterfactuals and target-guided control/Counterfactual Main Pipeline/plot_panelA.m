function plot_panelA(out_root, wave_id, varargin)
% PLOT_PANELA  (Panel A for Fig. X, v2: baseline line + translucent fills)
% 上图：f(t) —— 基线黑色粗实线；各场景半透明填充（“山峦图”）
% 下图：inc(t) —— 同样风格，并在右上角列出 %Avoided / Δinc_max / Δt_peak
%
% 依赖：smf_colorlib.m
%
% 用法：
%   plot_panelA('results_cf','Anhui-01', ...
%       'Scenarios', {'ADV3','NEG_STRONG','TGC_R2ZERO'}, ...
%       'SavePNG', true, 'FontSize', 11, 'LineWidth', 1.8);

% ---------- 参数 ----------
p = inputParser; p.FunctionName = mfilename;
p.addRequired('out_root', @(x)ischar(x)||isstring(x));
p.addRequired('wave_id',  @(x)ischar(x)||isstring(x));
p.addParameter('Scenarios', {'ADV3','NEG_STRONG','TGC_R2ZERO'}, ...
    @(x)iscellstr(x)||all(string(x)==string(x)));
p.addParameter('SavePNG',  false, @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('FontSize', 11,    @(x)isnumeric(x)&&isscalar(x)&&x>=6);
p.addParameter('LineWidth',1.8,   @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('OutName',  '',    @(x)ischar(x)||isstring(x));
p.addParameter('FillAlpha',0.28,  @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.parse(out_root, wave_id, varargin{:});
A = p.Results;

out_root = char(out_root);
wave_id  = char(wave_id);
fig_dir  = fullfile(out_root,'figures');
if ~isfolder(fig_dir), mkdir(fig_dir); end

% ---------- 读取结果 ----------
res_path = fullfile(out_root,'cf_results.mat');
if ~isfile(res_path), error('找不到结果文件：%s', res_path); end
S = load(res_path, 'Results', 'Scenarios');
Results = S.Results;

% ---------- 读取时序 ----------
ts_path = fullfile(out_root,'timeseries', sprintf('wave_%s.mat', sanitize(wave_id)));
if ~isfile(ts_path), error('找不到该波的时序：%s', ts_path); end
TSS = load(ts_path, 'TS'); TS = TSS.TS;

t    = TS.t(:);    % 注意：0..T-1
f0   = TS.f(:);
inc0 = TS.inc(:);
T    = numel(t);

% ---------- 场景列表 ----------
scns = string(A.Scenarios(:));
if isempty(scns)
    maskW = strcmp(string(Results.wave_id), string(wave_id));
    scns  = unique(string(Results.scenario(maskW)), 'stable');
end
Kscn = numel(scns);

% ---------- 配色 ----------
Cbase = [0 0 0];                                         % 基线黑
Cscn  = smf_colorlib('maina8','K', max(3, Kscn));        % 场景填充色
alpha = A.FillAlpha;

% ---------- 基线峰值与峰时 ----------
[inc0_max, idx0] = local_max_with_ties(inc0);
t0_peak = t(idx0);

% ---------- 统一结构体模板（修复“不同结构体下标赋值”报错） ----------
tmpl = struct('tag','','f',[],'inc',[],'pct',NaN,'dIncMax',NaN,'dTpeak',NaN);
Sinfo = repmat(tmpl, 0, 1);

% ---------- 收集各场景数据与指标 ----------
for i = 1:Kscn
    tag = upper(strtrim(scns(i)));
    f_field   = "f_"+tag;
    inc_field = "inc_"+tag;
    if ~isfield(TS, f_field) || ~isfield(TS, inc_field)
        warning('场景 %s：TS 中缺少 %s 或 %s，跳过。', tag, f_field, inc_field);
        continue;
    end
    f_s   = TS.(f_field)(:);
    inc_s = TS.(inc_field)(:);
    if numel(f_s)~=T || numel(inc_s)~=T
        warning('场景 %s：序列长度不符，跳过。', tag);
        continue;
    end

    % %Avoided：优先读 Results，否则即时计算
    pct = lookup_pct(Results, wave_id, tag);
    if isnan(pct)
        Cb = nansum(inc0); Cs = nansum(inc_s);
        pct = 100*(Cb - Cs)/max(Cb, eps);
    end

    % 峰值/峰时差
    [inc_s_max, idxs] = local_max_with_ties(inc_s);
    dIncMax = inc0_max - inc_s_max;
    dTpeak  = t(idxs) - t0_peak;

    % 填入（先扩容一格，再逐字段赋值，字段集合保持一致）
    Sinfo(end+1) = tmpl; %#ok<AGROW>
    Sinfo(end).tag   = char(tag);
    Sinfo(end).f     = f_s;
    Sinfo(end).inc   = inc_s;
    Sinfo(end).pct   = double(pct);
    Sinfo(end).dIncMax = double(dIncMax);
    Sinfo(end).dTpeak  = double(dTpeak);
end

if isempty(Sinfo), error('没有可用的场景时序可绘制。'); end

% ---------- 绘图 ----------
f = figure('Color','w'); 
tlo = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% ===== (a) f(t) =====
ax1 = nexttile(1); hold(ax1,'on'); box(ax1,'on');
plot(ax1, t, f0, '-', 'LineWidth', A.LineWidth+2.6, 'Color', [0 0 0], 'DisplayName','Baseline');

% crossing 与 t0 标注（若存在）
plot_meta_lines(ax1, TS);

% 叠加场景 f(t)
for i = 1:numel(Sinfo)
    ci = Cscn(i,:);
    plot(ax1, t, Sinfo(i).f, '-', 'LineWidth', A.LineWidth, 'Color', ci, ...
        'DisplayName', Sinfo(i).tag);
end
xlabel(ax1,'Day'); ylabel(ax1,'f(t)');
title(ax1, sprintf('%s — SMF edits', wave_id), 'Interpreter','none');
set(ax1,'FontSize',A.FontSize);
legend(ax1,'Location','best','Box','off'); grid(ax1,'on'); xlim(ax1,[min(t) max(t)]);


% ===== (b) inc(t) =====
ax2 = nexttile(2); hold(ax2,'on'); box(ax2,'on');

% 先画场景填充并收集句柄
hP = gobjects(numel(Sinfo),1);
for i = 1:numel(Sinfo)
    hP(i) = fill_between_zero(ax2, t, Sinfo(i).inc, Cscn(i,:), alpha);
end

% 再画基线并保存句柄
hBase = plot(ax2, t, inc0, '-', 'LineWidth', 3.0, 'Color', [0 0 0]);

% 加图例（基线 + 各场景）
legH = [hBase; hP(:)];
legT = ["Baseline"; string({Sinfo.tag})'];
legend(ax2, legH, legT, 'Location','best', 'Box','off', 'Interpreter','none');


xlabel(ax2,'Day'); ylabel(ax2,'inc(t) = \beta S I / N');
title(ax2, sprintf('%s — Incidence', wave_id), 'Interpreter','none');
set(ax2,'FontSize',A.FontSize); grid(ax2,'on'); xlim(ax2,[min(t) max(t)]);

% 右上角指标表
annotate_metrics(ax2, Sinfo, Cscn, A.FontSize);

% ---------- 保存 ----------
if A.SavePNG
    if A.OutName==""
        out_png = fullfile(fig_dir, sprintf('panelA_%s.png', sanitize(wave_id)));
    else
        out_png = fullfile(fig_dir, char(A.OutName));
    end
    exportgraphics(tlo, out_png, 'Resolution', 300);
    fprintf('[plot_panelA] Saved → %s\n', out_png);
end
end % ===== 主函数 =====


% ======== 工具函数 ========
function s = sanitize(s), s = regexprep(char(s), '[^a-zA-Z0-9-_]', '_'); end

function [mx, idx] = local_max_with_ties(y)
y = y(:);
[~, idx] = max(y); 
mx = y(idx);
end

function pct = lookup_pct(Results, wave_id, tag)
pct = NaN;
if isempty(Results), return; end
mask = strcmpi(string(Results.wave_id), string(wave_id)) & ...
       strcmpi(string(Results.scenario), string(tag));
if any(mask)
    v = Results.avoided_pct(mask);
    if ~isempty(v), pct = double(v(1)); end
end
end

function plot_meta_lines(ax, TS)
if isfield(TS,'meta') && isstruct(TS.meta)
    if isfield(TS.meta,'t_cross') && ~isempty(TS.meta.t_cross) && ~isnan(TS.meta.t_cross)
        xline(ax, TS.meta.t_cross - 1, ':', 'Color', [0.5 0.5 0.5], ...
              'LineWidth', 2.4, 'HandleVisibility','off');
    end
    if isfield(TS.meta,'t0') && ~isempty(TS.meta.t0) && ~isnan(TS.meta.t0)
        xline(ax, TS.meta.t0 - 1, '--', 'Color', [0.85 0.2 0.2], ...
              'LineWidth', 2.4, 'HandleVisibility','off');
    end
end
end


function h = fill_between_zero(ax, x, y, col, alpha)
x = x(:); y = y(:);
X = [x; flipud(x)];
Y = [zeros(size(y)); flipud(y)];
h = patch('Parent',ax, 'XData',X, 'YData',Y, ...
          'FaceColor',col, 'EdgeColor',[0 0 0], 'FaceAlpha',alpha);
end




function annotate_metrics(ax, Sinfo, Cscn, fs)
xl = xlim(ax); yl = ylim(ax);
x0 = xl(1) + 0.70*(xl(2)-xl(1));
y0 = yl(1) + 0.98*(yl(2)-yl(1));
dy = 0.06*(yl(2)-yl(1));  % 行距
text(ax, x0, y0, " %Avoided   \Delta inc_{max}   \Delta t_{peak}", ...
    'HorizontalAlignment','left','VerticalAlignment','top', ...
    'FontSize', fs, 'Color', [0.3 0.3 0.3], 'Interpreter','tex');

for i = 1:numel(Sinfo)
    yi = y0 - i*dy;
    s  = sprintf(' %6.1f%%    %8.3g     %+6.0f', Sinfo(i).pct, Sinfo(i).dIncMax, Sinfo(i).dTpeak);
    text(ax, x0, yi, [Sinfo(i).tag, ' |', s], ...
        'HorizontalAlignment','left','VerticalAlignment','top', ...
        'FontSize', fs, 'Color', Cscn(i,:), 'Interpreter','tex');
end
end

