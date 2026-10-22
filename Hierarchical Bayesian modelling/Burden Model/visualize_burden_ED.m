function visualize_burden_ED(varargin)
% VISUALIZE_BURDEN_ED (robust)
% 一键生成扩展图 ED1–ED3（省级/变异株随机截距、PPC 覆盖率）
%
% 用法：
%   visualize_burden_ED('results_root','results_hbm','dpi',400)

p=inputParser;
addParameter(p,'results_root','results_hbm',@(s)ischar(s)||isstring(s));
addParameter(p,'dpi',400,@(x)isnumeric(x)&&x>0);
addParameter(p,'font','Arial',@(s)ischar(s)||isstring(s));
parse(p,varargin{:}); args=p.Results;

root = char(args.results_root);
burden_dir = fullfile(root,'burden');
if ~isfolder(burden_dir), error('缺少目录：%s', burden_dir); end
out_dir = fullfile(root,'figs_burden'); if ~isfolder(out_dir), mkdir(out_dir); end

% ---------- 识别两个目标子目录（兼容 *_obs 命名） ----------
[tgtI, tgtC] = pick_targets(burden_dir);
fprintf('[ED] 使用子目录：I_max=%s | C_total=%s\n', tgtI, tgtC);

% 统一配色
COL.gray = [0.15 0.15 0.15];
COL.grid = [0.90 0.90 0.90];
COL.variants = [0.82 0.33 0.30; 0.95 0.66 0.12; 0.12 0.47 0.71];


[prov_names_map, var_names_map] = build_label_maps(root);

%% ============ ED1：省级随机截距毛毛虫（两面板） ============
f1 = figure('Color','w','Units','pixels','Position',[80 80 620 900]);
tiledlayout(f1,2,1,'Padding','compact','TileSpacing','compact');

% I_max
nexttile; 
plot_caterpillar(fullfile(burden_dir,tgtI,'ranef_province.csv'), ...
    'Province random intercepts: I_{max}', COL.gray, args.font, prov_names_map);

nexttile; 
plot_caterpillar(fullfile(burden_dir,tgtC,'ranef_province.csv'), ...
    'Province random intercepts: C_{total}', COL.gray, args.font, prov_names_map);

exportgraphics(f1, fullfile(out_dir,'ED1_province_ranef.png'), 'Resolution', args.dpi);
print(f1, fullfile(out_dir,'ED1_province_ranef.pdf'), '-dpdf','-painters');

%% ============ ED2：变异株随机截距比较（两面板） ============
f2 = figure('Color','w','Units','pixels','Position',[80 80 900 380]);
tiledlayout(f2,1,2,'Padding','compact','TileSpacing','compact');

% I_max
nexttile;
plot_variant_re(fullfile(burden_dir,tgtI,'ranef_variant.csv'), ...
    'Variant random intercepts: I_{max}', COL.variants, COL.gray, args.font, var_names_map);

nexttile;
plot_variant_re(fullfile(burden_dir,tgtC,'ranef_variant.csv'), ...
    'Variant random intercepts: C_{total}', COL.variants, COL.gray, args.font, var_names_map);

exportgraphics(f2, fullfile(out_dir,'ED2_variant_ranef.png'), 'Resolution', args.dpi);
print(f2, fullfile(out_dir,'ED2_variant_ranef.pdf'), '-dpdf','-painters');

%% ============ ED3：PPC 覆盖率摘要（两面板） ============
f3 = figure('Color','w','Units','pixels','Position',[80 80 900 380]);
tiledlayout(f3,1,2,'Padding','compact','TileSpacing','compact');

% I_max
nexttile;
plot_ppc_coverage_bar(fullfile(burden_dir,tgtI,'ppc.csv'), ...
    'PPC coverage: I_{max}', args.font);

% C_total
nexttile;
plot_ppc_coverage_bar(fullfile(burden_dir,tgtC,'ppc.csv'), ...
    'PPC coverage: C_{total}', args.font);

exportgraphics(f3, fullfile(out_dir,'ED3_ppc_coverage.png'), 'Resolution', args.dpi);
print(f3, fullfile(out_dir,'ED3_ppc_coverage.pdf'), '-dpdf','-painters');

fprintf('[ED] 完成。输出目录：%s\n', out_dir);
end

% ======================== 工具函数 ========================

function [tgtI, tgtC] = pick_targets(burden_dir)
d=dir(burden_dir); d=d([d.isdir]); names=setdiff({d.name},{'.','..'});
cands = {};
for i=1:numel(names)
    if isfile(fullfile(burden_dir,names{i},'summary.csv'))
        cands{end+1}=names{i}; %#ok<AGROW>
    end
end
if numel(cands)<2
    error('未找到足够的目标子目录（包含 summary.csv）。发现：%s', strjoin(cands,','));
end
% 尝试匹配
i_idx = find(contains(lower(cands),'i'),1,'first'); 
c_idx = find(contains(lower(cands),'c'),1,'first'); 
if isempty(i_idx) || isempty(c_idx)
    % 兜底：前两个
    tgtI = cands{1}; tgtC = cands{2};
else
    tgtI = cands{i_idx}; tgtC = cands{c_idx};
end
end

function T = read_csv_robust(csv_path)
try
    opts = detectImportOptions(csv_path);
    opts.VariableNamingRule = 'preserve';
    T = readtable(csv_path, opts);
catch
    T = readtable(csv_path, 'VariableNamingRule','preserve');
end
end


function v = col_vec(T, candidates)
% 在表 T 中不区分大小写地匹配候选列名（去掉空格/下划线）
vn = string(T.Properties.VariableNames);
vn_clean = lower(erase(erase(vn,"_")," "));
cand_clean = lower(erase(erase(string(candidates),"_")," "));
idx = [];
for c = cand_clean(:)'
    k = find(vn_clean==c, 1, 'first');
    if ~isempty(k), idx = k; break; end
end
if isempty(idx)
    error('列缺失。需要其中之一：%s | 实际列名：%s', ...
        strjoin(string(candidates),', '), strjoin(vn,', '));
end
v = T{:, idx};
end







function plot_variant_re(csv_path, ttl, cols, col_line, fontname, var_names_map)
if ~isfile(csv_path)
    text(0.5,0.5,sprintf('Missing: %s',csv_path),'HorizontalAlignment','center'); axis off; return;
end
R = read_csv_robust(csv_path);

mu  = try_cols(R, {'mean','Mean','MEAN','mu','avg','effect','Effect','EFFECT'});
ciL = try_cols(R, {'ci2_5','CI2_5','ci.2.5','lower','lcl','lci'}, true);
ciU = try_cols(R, {'ci97_5','CI97_5','ci97.5','upper','ucl','uci'}, true);

% 标签：CSV 的 variant；否则 variant_id → 映射；再否则兜底
if any(strcmpi(R.Properties.VariableNames,'variant'))
    vlab = string(R{:, find(strcmpi(R.Properties.VariableNames,'variant'),1,'first')});
elseif any(strcmpi(R.Properties.VariableNames,'variant_id')) && ~isempty(var_names_map)
    vid = try_cols(R, {'variant_id','var_id','var'});
    vid = max(1, min(numel(var_names_map), round(vid(:))));
    vlab = var_names_map(vid);
else
    vlab = string(arrayfun(@(k)sprintf('Variant %d',k), 1:height(R), 'UniformOutput',false));
end
vlab = string(vlab(:));

mask = isfinite(mu);
if ~isempty(ciL), mask=mask & isfinite(ciL); end
if ~isempty(ciU), mask=mask & isfinite(ciU); end
mu=mu(mask); vlab=vlab(mask);
if ~isempty(ciL), ciL=ciL(mask); ciU=ciU(mask); end
if isempty(mu), text(0.5,0.5,'No data','HorizontalAlignment','center'); axis off; return; end

hold on; box off; grid on;
for i=1:numel(mu)
    if ~isempty(ciL)
        plot([ciL(i) ciU(i)], [i i], '-', 'Color', col_line, 'LineWidth',1.2);
    end
    c = cols( 1+mod(i-1,size(cols,1)) , : );
    plot(mu(i), i, 'o', 'MarkerFaceColor',c,'MarkerEdgeColor','k','MarkerSize',6);
end
xline(0,'-','Color',[0.6 0.6 0.6],'LineWidth',1);

yticks(1:numel(mu)); yticklabels(vlab);
set(gca,'YDir','reverse','TickDir','out','LineWidth',0.8,'FontName',fontname,'FontSize',10);
xlabel('Random intercept','FontName',fontname,'FontSize',11);
title(ttl,'FontWeight','normal','FontName',fontname,'FontSize',12);


ax=gca; ax.GridColor=[0.9 0.9 0.9]; ax.GridAlpha=0.85;
if ~isempty(ciL), lo=min(ciL); hi=max(ciU); else, lo=min(mu); hi=max(mu); end
if ~(hi>lo), lo=lo-1; hi=hi+1; end
pad=0.05*(hi-lo); xlim([lo-pad, hi+pad]);
end









function plot_caterpillar(csv_path, ttl, col_line, fontname, prov_names_map)
if ~isfile(csv_path)
    text(0.5,0.5,sprintf('Missing: %s',csv_path),'HorizontalAlignment','center'); axis off; return;
end
R = read_csv_robust(csv_path);

mu  = try_cols(R, {'mean','Mean','MEAN','mu','avg','effect','Effect','EFFECT'});
ciL = try_cols(R, {'ci2_5','CI2_5','ci.2.5','lower','lcl','lci'}, true);
ciU = try_cols(R, {'ci97_5','CI97_5','ci97.5','upper','ucl','uci'}, true);

% 省份标签：优先 CSV 自带 name/province，其次 province_id → 映射
labs = [];
if any(strcmpi(R.Properties.VariableNames,'name'))
    labs = string(R{:, find(strcmpi(R.Properties.VariableNames,'name'),1,'first')});
elseif any(strcmpi(R.Properties.VariableNames,'province'))
    labs = string(R{:, find(strcmpi(R.Properties.VariableNames,'province'),1,'first')});
end
if isempty(labs) && any(strcmpi(R.Properties.VariableNames,'province_id')) && ~isempty(prov_names_map)
    pid = try_cols(R, {'province_id','prov_id','prov'});
    pid = max(1, min(numel(prov_names_map), round(pid(:))));
    labs = prov_names_map(pid);
end
if isempty(labs)
    labs = arrayfun(@(k)sprintf('Prov %d',k), 1:height(R), 'UniformOutput',false);
end
labs = string(labs(:));

mask = isfinite(mu);
if ~isempty(ciL), mask = mask & isfinite(ciL); end
if ~isempty(ciU), mask = mask & isfinite(ciU); end
mu=mu(mask); labs=labs(mask);
if ~isempty(ciL), ciL=ciL(mask); ciU=ciU(mask); end
if isempty(mu), text(0.5,0.5,'No data','HorizontalAlignment','center'); axis off; return; end

[~,ord] = sort(mu);
mu=mu(ord); labs=labs(ord);
if ~isempty(ciL), ciL=ciL(ord); ciU=ciU(ord); end

y = 1:numel(mu); hold on; box off; grid on;
if ~isempty(ciL)
    plot([ciL ciU].', [y;y], '-', 'Color', col_line, 'LineWidth',1.2);
end
plot(mu, y, 'o', 'MarkerFaceColor',[0.1 0.1 0.1], 'MarkerEdgeColor','w','MarkerSize',4);
xline(0,'-','Color',[0.6 0.6 0.6],'LineWidth',1);

set(gca,'YTick',y,'YTickLabel',labs,'YDir','reverse',...
    'TickDir','out','LineWidth',0.8,'FontName',fontname,'FontSize',10);
xlabel('Random intercept','FontName',fontname,'FontSize',11);
title(ttl,'FontWeight','normal','FontName',fontname,'FontSize',12);


ax=gca; ax.GridColor=[0.9 0.9 0.9]; ax.GridAlpha=0.85;
if ~isempty(ciL), lo=min(ciL); hi=max(ciU); else, lo=min(mu); hi=max(mu); end
if ~(hi>lo), lo=lo-1; hi=hi+1; end
pad=0.05*(hi-lo); xlim([lo-pad, hi+pad]);
end









function plot_ppc_coverage_bar(ppc_csv, ttl, fontname)
if ~isfile(ppc_csv)
    text(0.5,0.5,sprintf('Missing: %s',ppc_csv),'HorizontalAlignment','center'); axis off; return;
end
T = read_csv_robust(ppc_csv);

has10  = all(ismember({'yhat_p10','yhat_p90'}, T.Properties.VariableNames));
has025 = all(ismember({'yhat_p2_5','yhat_p97_5'}, T.Properties.VariableNames));

y = col_vec(T, {'y_obs','yobs','obs'});
ok = isfinite(y);
if has10,  ok = ok & isfinite(col_vec(T,{'yhat_p10'})) & isfinite(col_vec(T,{'yhat_p90'})); end
if has025, ok = ok & isfinite(col_vec(T,{'yhat_p2_5'})) & isfinite(col_vec(T,{'yhat_p97_5'})); end
T = T(ok,:);

cov10 = NaN; cov025 = NaN;
if has10
    y10 = col_vec(T,{'yhat_p10'}); y90 = col_vec(T,{'yhat_p90'});
    cov10  = mean( (col_vec(T,{'y_obs'})>=y10) & (col_vec(T,{'y_obs'})<=y90) );
end
if has025
    y025 = col_vec(T,{'yhat_p2_5'}); y975 = col_vec(T,{'yhat_p97_5'});
    cov025 = mean( (col_vec(T,{'y_obs'})>=y025) & (col_vec(T,{'y_obs'})<=y975) );
end

vals = [cov10, cov025];
labs = {'10–90%','2.5–97.5%'};
clr  = [0.25 0.47 0.78; 0.30 0.70 0.52];

hold on; box off; grid on;
for i=1:numel(vals)
    if ~isnan(vals(i))
        bar(i, vals(i), 0.55, 'FaceColor', clr(i,:), 'EdgeColor','k','LineWidth',0.8);
        text(i, min(0.98, vals(i)+0.025), sprintf('%.3f', vals(i)), 'HorizontalAlignment','center', ...
            'FontName',fontname,'FontSize',10);
    end
end
set(gca,'XTick',1:numel(vals),'XTickLabel',labs,'YLim',[0 1],...
    'TickDir','out','LineWidth',0.8,'FontName',fontname,'FontSize',10);
ylabel('Coverage','FontName',fontname,'FontSize',11);
title(ttl,'FontWeight','normal','FontName',fontname,'FontSize',12);

ax=gca; ax.GridColor=[0.9 0.9 0.9]; ax.GridAlpha=0.85;
end


function v = try_cols(T, candidates, allowEmpty)
if nargin<3, allowEmpty=false; end
vn = string(T.Properties.VariableNames);
vn_clean = lower(erase(erase(vn,"_")," "));
cand_clean = lower(erase(erase(string(candidates),"_")," "));
idx = [];
for c = cand_clean(:)'
    k = find(vn_clean==c, 1, 'first');
    if ~isempty(k), idx = k; break; end
end
if isempty(idx)
    if allowEmpty
        v = [];
    else
        error('列缺失。需要其中之一：%s | 实际列名：%s', ...
            strjoin(string(candidates),', '), strjoin(vn,', '));
    end
else
    v = T{:, idx};
end
end

function [prov_names, var_names] = build_label_maps(results_root)
% 返回：prov_names(Px1 string), var_names(Vx1 string)
% 逻辑：读 cache 获取 id；从 features_aug.mat（或 features.mat）找名字列；
%      对每个 id 取众数文本；拿不到就退回 "Prov #"/"Variant #"

% 读 cache
C = load(fullfile(results_root,'cache','hbm_cache.mat')); cache = C.cache;
P = cache.P; V = cache.V;

% 找 features 文件
feat = '';
cand = { fullfile(results_root,'..','features_aug.mat'), ...
         fullfile(results_root,'..','features.mat'), ...
         'features_aug.mat', 'features.mat' };
for i=1:numel(cand)
    if isfile(cand{i}), feat = cand{i}; break; end
end
prov_names = "Prov " + string((1:P)');
var_names  = "Variant " + string((1:V)');

if isempty(feat), return; end

S = load(feat);
T = [];
if isfield(S,'Features') && istable(S.Features), T = S.Features; end
if isempty(T), return; end

% —— 省名列候选
prov_col = pick_text_col(T, {'province','Province','prov_name','prov','ProvinceName','region'});
% —— 毒株列候选
var_col  = pick_text_col(T, {'variant','Variant','lineage','VOC','voc'});

if ~isempty(prov_col)
    prov_text = string(T.(prov_col));
    prov_text(prov_text=="") = missing;   % 空变 missing 才能正确众数
    prov_names = strings(P,1);
    for id=1:P
        ix = (cache.prov_id==id);
        prov_names(id) = safe_mode(prov_text(ix), "Prov " + string(id));
    end
end

if ~isempty(var_col)
    var_text = string(T.(var_col));
    var_text(var_text=="") = missing;
    var_names = strings(V,1);
    for id=1:V
        ix = (cache.var_id==id);
        var_names(id) = safe_mode(var_text(ix), "Variant " + string(id));
    end
end
end

function colname = pick_text_col(T, candidates)
colname = "";
v = string(T.Properties.VariableNames);
% 先按候选名匹配
for c = string(candidates)
    k = find(strcmpi(v,c),1,'first');
    if ~isempty(k), colname = v(k); return; end
end
% 再退回寻找首个文本/分类列
for k=1:width(T)
    x = T{:,k};
    if iscellstr(x) || isstring(x) || iscategorical(x)
        colname = v(k); return;
    end
end
end

function s = safe_mode(vec, fallback)
% 求众数（忽略 missing）；全 missing 则用 fallback
if iscategorical(vec), vec = string(vec); end
vec = vec(~ismissing(vec));
if isempty(vec), s = fallback; return; end
u = unique(vec);
counts = zeros(numel(u),1);
for i=1:numel(u), counts(i)=sum(vec==u(i)); end
[~,m] = max(counts);
s = u(m);
end
