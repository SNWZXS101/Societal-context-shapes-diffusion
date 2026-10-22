function ProvScores = residualize_value_added(WaveScores, varargin)
% RESIDUALIZE_VALUE_ADDED
% 省级价值增量（VA）与综合 SRS 的稳健估计：
% - RRI_VA    = 以 RRI_wave 为因变量的省随机截距（或回退 WLS）→ 0–100
% - MEI_VA    = 以 MEI_wave 为因变量的省随机截距（或回退 WLS）→ 0–100
% - SRS       = w_RRI * RRI_VA + w_MEI * MEI_VA   （默认 0.4 / 0.6；可调）
% - 可选：引入协变量 X（固定效应），以及变异株固定效应 variant
% - 采用权重列 weights_col（默认 'w_tail'；找不到时报错）
%
% 输出：
%   <outdir>/map_ready_SRS.csv（无时间戳，便于下游脚本固定读取）
%   <outdir>/map_ready_SRS_RRIMEI_YYYYMMDD_HHMM.csv（带时间戳的完整表）
%
% 例：
% ProvScores = residualize_value_added(WaveScores, ...
%     'use_variant',true, ...
%     'xvars',{'urban_rate','log_income_pc','log_transit_vol'}, ...
%     'features_mat','features_aug.mat', ...
%     'weights_col','w_tail', ...
%     'w_RRI',0.4,'w_MEI',0.6, ...
%     'outdir','scores','verbose',true);

% ---------- 参数 ----------
p = inputParser; p.FunctionName = mfilename;
addRequired(p,'WaveScores',@(t) istable(t));
addParameter(p,'scenario','',@(s) ischar(s)||isstring(s));   % 容错：若 WaveScores 有 scenario 列
addParameter(p,'use_variant',true,@(x) islogical(x)||ismember(x,[0 1]));
addParameter(p,'weights_col','w_tail',@(s) ischar(s)||isstring(s));
addParameter(p,'outdir','scores',@(s) ischar(s)||isstring(s));
addParameter(p,'verbose',true,@(x) islogical(x)||ismember(x,[0 1]));
% 新增：协变量接口（可选）
addParameter(p,'xvars',{},@(c) iscellstr(c) || all(arrayfun(@(z) isstring(z)||ischar(z), c)));
addParameter(p,'features_mat','',@(s) ischar(s)||isstring(s));
% 新增：SRS 权重
addParameter(p,'w_RRI',0.4,@(x) isnumeric(x)&&isscalar(x));
addParameter(p,'w_MEI',0.6,@(x) isnumeric(x)&&isscalar(x));
% 新增：是否落带时间戳的完整表
addParameter(p,'timestamp',true,@(x) islogical(x)||ismember(x,[0 1]));
parse(p,WaveScores,varargin{:});
A = p.Results;
if ~isfolder(A.outdir), mkdir(A.outdir); end

% ---- 取波级数据（如存在 scenario 列可筛选）----
T = WaveScores;
if ~isempty(A.scenario) && ismember('scenario', T.Properties.VariableNames)
    T = T(strcmpi(string(T.scenario), string(A.scenario)), :);
end

% ---- 必要列检查 ----
needBase = {'province','wave_id'};
assert(all(ismember(needBase, T.Properties.VariableNames)), 'WaveScores 缺少必要列 province/wave_id');

% ---- 权重列 ----
wcol = string(A.weights_col);
if ~ismember(wcol, T.Properties.VariableNames)
    error('未找到权重列：%s', wcol);
end

% ---- 变异株列（可选）----
hasVarCol = A.use_variant && ismember('variant', T.Properties.VariableNames);
if hasVarCol, T.variant = categorical(string(T.variant)); end

% ---- 统一类型/基础清洗 ----
T.province = categorical(string(T.province));
T = T(isfinite(T.(wcol)) & T.(wcol) > 0, :);

% 每省≥2波（保证随机效应可估）
G = groupsummary(T,'province');
drop1 = string(G.province(G.GroupCount<=1));
if ~isempty(drop1) && A.verbose
    fprintf('[residualize_value_added] 仅1波省剔除：%s\n', strjoin(cellstr(drop1),', '));
end
if ~isempty(drop1), T = T(~ismember(string(T.province), drop1), :); end

% ---- 协变量并入（可选）----
X_terms = "";
X_tbl = table();
if ~isempty(A.xvars)
    if isempty(A.features_mat) || ~isfile(A.features_mat)
        error('指定了 xvars 但未提供可用的 features_mat 文件。');
    end
    S = load(A.features_mat);
    Tfeat = [];
    for nm = {'Features','Tfeat','features','T'}
        if isfield(S, nm{1}), Tfeat = S.(nm{1}); break; end
    end
    if ~istable(Tfeat)
        error('features_mat 中未找到表（Features/Tfeat/features/T）。');
    end
    % 仅保留需要的键与协变量列
    keys = {'province','wave_id'};
    miss = setdiff(A.xvars, Tfeat.Properties.VariableNames);
    if ~isempty(miss)
        error('features_mat 缺少协变量列：%s', strjoin(miss, ', '));
    end
    X_tbl = Tfeat(:, [intersect(keys, Tfeat.Properties.VariableNames), A.xvars]);
    % 统一类型
    if ~isstring(X_tbl.province), X_tbl.province = string(X_tbl.province); end
    if ~isstring(X_tbl.wave_id),  X_tbl.wave_id  = string(X_tbl.wave_id);  end
    T = outerjoin(T, X_tbl, 'Keys',{'province','wave_id'}, 'MergeKeys',true);
    % 标准化各 X（z-score；忽略缺失）
    for j = 1:numel(A.xvars)
        v = A.xvars{j};
        if ~ismember(v, T.Properties.VariableNames), error('缺少并入后的协变量：%s', v); end
        x = double(T.(v));
        mu = mean(x,'omitnan'); sg = std(x,'omitnan');
        if ~isfinite(sg) || sg==0, sg = 1; end
        T.(v) = (x - mu) ./ sg;
    end
    % 生成公式片段
    X_terms = strjoin(string(A.xvars), ' + ');
end

% ---------- 拟合器：给定 target（RRI_wave / MEI_wave） ----------
    function OUT = fit_target(targetName)
        OUT = struct('VA',table(),'choice','', 'info_with','', 'info_novar','');
        if ~ismember(targetName, T.Properties.VariableNames)
            OUT.VA = table(string([]),[],[], 'VariableNames',{'province','VA_raw','n_obs'});
            OUT.choice = 'missing';
            return;
        end
        S = T(isfinite(T.(targetName)), :);
        if isempty(S)
            OUT.VA = table(string([]),[],[], 'VariableNames',{'province','VA_raw','n_obs'});
            OUT.choice = 'empty';
            return;
        end
        wts = double(S.(wcol));

        % 子函数：一套模型（含/不含变异株），优先 LME，失败则 WLS 回退
                function [VA_tbl, deg, info] = fit_one(spec_hasVariant)
            % ------- 始终先用 WLS 稳健回归 -------
            info = struct('path','WLS','hasVariant',spec_hasVariant,'msg','');
            y = S.(targetName);

            % 1) 组装设计矩阵 X（截距 + Xvars + variant dummies）
            X = ones(height(S),1);  % (Intercept)
            if ~isempty(A.xvars)
                for jj=1:numel(A.xvars)
                    xv = double(S.(A.xvars{jj}));
                    xv(~isfinite(xv)) = 0; % 缺失→0（等于均值填补后的标准化）
                    X = [X, xv]; %#ok<AGROW>
                end
            end
            if spec_hasVariant && hasVarCol
                [~,~,vlabs] = unique(S.variant);
                K = max(vlabs);
                if K>1
                    Xv = zeros(numel(y), K-1);
                    for j=2:K, Xv(vlabs==j, j-1)=1; end
                    X = [X, Xv]; %#ok<AGROW>
                end
            end

            % 2) 加权最小二乘（lscov）
            wts = double(S.(wcol));
            [beta, ~] = lscov(X, y, wts);

            % 3) 剥离固定效应，取省内加权均值 = 省“价值增量”
            if size(X,2) > 1
                y_adj = y - X(:,2:end) * beta(2:end);
            else
                y_adj = y;
            end
            P = string(S.province); up = unique(P);
            VA_raw = nan(numel(up),1); nobs = zeros(numel(up),1);
            for i=1:numel(up)
                m = (P==up(i)); ww=wts(m); yy=y_adj(m);
                VA_raw(i) = sum(ww.*yy)/sum(ww); nobs(i) = sum(m);
            end
            VA_tbl = table(up, VA_raw, nobs, 'VariableNames', {'province','VA_raw','n_obs'});

            % 4) 退化检测
            v = VA_tbl.VA_raw;
            deg = ~any(isfinite(v)) || range(v(~isnan(v))) < 1e-8 ...
                  || numel(unique(round(v(~isnan(v)),12)))<=1;

            % 5) 如确有必要，再尝试 LME 获取 RE（失败也不报错，只做静默回退）
            %    —— 大多数情况下 WLS 已足够且更稳健
            if deg==false
                return
            end
            try
                rhs_terms = ["1"];
                if ~isempty(A.xvars), rhs_terms = [rhs_terms, string(A.xvars)]; end %#ok<AGROW>
                if spec_hasVariant && hasVarCol, rhs_terms = [rhs_terms, "variant"]; end %#ok<AGROW>
                rhs = strjoin(rhs_terms, " + ");
                formula = sprintf('%s ~ %s + (1|province)', targetName, rhs);
                lme = fitlme(S, formula, 'FitMethod','REML','DummyVarCoding','effects','Weights',wts);

                RE = randomEffects(lme);  % R2024b 通常返回 table
                if isa(RE,'dataset'), RE = dataset2table(RE); end
                if isstruct(RE), RE = struct2table(RE); end

                % 宽松匹配列名（不同版本可能差异）
                cols = lower(string(RE.Properties.VariableNames));
                kName   = find(cols=="name",1);
                kGroup  = find(cols=="group" | cols=="groupingvariable" | cols=="groupname",1);
                kLevel  = find(cols=="level",1);
                kEst    = find(cols=="estimate" | cols=="value",1);
                if ~isempty(kEst) && ~isempty(kLevel)
                    namecol  = RE{:, kName};   % 可能不存在
                    if isempty(namecol), namecol = repmat("(Intercept)", height(RE),1); end
                    isProv = true(size(RE,1),1);
                    if ~isempty(kGroup)
                        isProv = isProv & strcmpi(string(RE{:,kGroup}),"province");
                    end
                    isInt  = isProv & strcmpi(string(namecol),"(Intercept)");
                    lvl    = string(RE{isInt, kLevel});
                    est    = double(RE{isInt, kEst});
                    if ~isempty(lvl)
                        VA_tbl = table(lvl, est, 'VariableNames', {'province','VA_raw'});
                        % 补 n_obs
                        Nobs = groupsummary(S,'province');
                        nmap = containers.Map(cellstr(string(Nobs.province)), num2cell(Nobs.GroupCount));
                        VA_tbl.n_obs = zeros(height(VA_tbl),1);
                        for ii=1:height(VA_tbl)
                            k = char(VA_tbl.province(ii)); if isKey(nmap,k), VA_tbl.n_obs(ii) = nmap(k); end
                        end
                        deg = false; info.path = 'LME';
                    end
                end
            catch
                % 保持 WLS 结果；不改变 info，也不抛错
            end
        end


        [VA_var, deg_var, info_var] = fit_one(true);
        [VA_nov, deg_nov, info_nov] = fit_one(false);

        if ~deg_var
            OUT.VA = VA_var; OUT.choice = 'with_variant';
        else
            OUT.VA = VA_nov; OUT.choice = 'no_variant';
        end
        OUT.info_with  = human_info(info_var);
        OUT.info_novar = human_info(info_nov);
    end

% --- 分别拟合两口径 ---
OUT_RRI = fit_target('RRI_wave');
OUT_MEI = fit_target('MEI_wave');

% --- 重命名后合并 ---
R = OUT_RRI.VA; E = OUT_MEI.VA;
if ~isempty(R), R.province = string(R.province); R = renamevars(R, {'VA_raw','n_obs'}, {'VA_raw_RRI','n_obs_RRI'}); end
if ~isempty(E), E.province = string(E.province); E = renamevars(E, {'VA_raw','n_obs'}, {'VA_raw_MEI','n_obs_MEI'}); end
if ~isempty(R) && ~isempty(E)
    ProvScores = outerjoin(R, E, 'Keys','province','MergeKeys',true,'Type','full');
elseif ~isempty(R)
    ProvScores = R; ProvScores.VA_raw_MEI = NaN(height(ProvScores),1); ProvScores.n_obs_MEI = NaN(height(ProvScores),1);
elseif ~isempty(E)
    ProvScores = E; ProvScores.VA_raw_RRI = NaN(height(ProvScores),1); ProvScores.n_obs_RRI = NaN(height(ProvScores),1);
else
    ProvScores = table(string([]), [], [], [], [], 'VariableNames', ...
        {'province','VA_raw_RRI','n_obs_RRI','VA_raw_MEI','n_obs_MEI'});
end

% --- 0–100 映射：RRI_VA / MEI_VA，并合成 SRS ---
ProvScores.RRI_VA = map01_rank(ProvScores.VA_raw_RRI);
ProvScores.MEI_VA = map01_rank(ProvScores.VA_raw_MEI);

% 若有全缺省或常数退化，map01_rank 会回到 50 分常数
w_sum = A.w_RRI + A.w_MEI; if w_sum==0, w_sum = 1; end
ProvScores.SRS = (A.w_RRI * ProvScores.RRI_VA + A.w_MEI * ProvScores.MEI_VA) / w_sum;

% 模型选择信息
ProvScores.model_choice_RRI = repmat(string(OUT_RRI.choice), height(ProvScores),1);
ProvScores.model_choice_MEI = repmat(string(OUT_MEI.choice), height(ProvScores),1);

% 排序与输出（按 SRS 降序）
[~,rk] = sort(ProvScores.SRS,'descend','MissingPlacement','last');
ProvScores = ProvScores(rk,:);
ProvScores.rank = (1:height(ProvScores))';

% --- 落盘：一个固定名 + 一个带时间戳完整表 ---
outcsv_fixed = fullfile(A.outdir,'map_ready_SRS.csv');
writetable(ProvScores(:, {'province','RRI_VA','MEI_VA','SRS','rank'}), outcsv_fixed);

tag = datestr(now,'yyyymmdd_HHMM');
if A.timestamp
    suffix = "_" + string(tag);
else
    suffix = "";
end
outcsv_full = fullfile(A.outdir, sprintf('map_ready_SRS_RRIMEI%s.csv', suffix));

writetable(ProvScores, outcsv_full);

if A.verbose
    fprintf('[residualize_value_added] 省数=%d | 选择：RRI=%s / MEI=%s → %s\n', ...
        height(ProvScores), OUT_RRI.choice, OUT_MEI.choice, outcsv_fixed);
    fprintf('  RRI  withVar:%s | noVar:%s\n', OUT_RRI.info_with, OUT_RRI.info_novar);
    fprintf('  MEI  withVar:%s | noVar:%s\n', OUT_MEI.info_with, OUT_MEI.info_novar);
    fprintf('  SRS 权重：w_RRI=%.2f, w_MEI=%.2f | 完整输出：%s\n', A.w_RRI, A.w_MEI, outcsv_full);
end
end

% ---------- helpers ----------
function s = human_info(z)
if isstruct(z)
    s = sprintf('%s%s', z.path, tern(z.hasVariant,'+variant',''));
    if isfield(z,'msg') && ~isempty(z.msg), s = [s, ' (', z.msg, ')']; end
else
    s = "";
end
end
function y = tern(cond,a,b), if cond, y=a; else, y=b; end, end
function s = map01_rank(v)
v = double(v(:));
if isempty(v)
    s = v; return;
end
if all(~isfinite(v)) || numel(unique(round(v(~isnan(v)),12)))<=1
    s = 50*ones(size(v));
else
    r = tiedrank(v); s = 100*(r - 0.5)/numel(r);
end
end
