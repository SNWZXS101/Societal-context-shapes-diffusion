function fit_hbm_main_m0_matlab(varargin)
% 分层 Student-t 回归（纯 MATLAB 版本）
% 支持：二次项、hinge样条、异方差、随机斜率、可调先验强度
% ---- 调用示例 ----
% fit_hbm_main_m0_matlab('results_root','results_hbm', ...
%   'num_warmup',4000,'num_samples',4000,'chains',4,'seed',2025, ...
%   'beta_sd',2,'alpha_sd',3,'ranef_scale',3,'learn_nu',false, ...
%   'nu_step',0.10,'init_nu',30,'init_jitter',0.3, ...
%   'rand_slope_prov_idx',[7 12], 'heteroskedastic_idx',[8]);

% ================== 解析参数 ==================
p = inputParser;
addParameter(p,'results_root','results_hbm',@ischar);
addParameter(p,'num_warmup',1000,@(x)isnumeric(x)&&x>=10);
addParameter(p,'num_samples',1000,@(x)isnumeric(x)&&x>=10);
addParameter(p,'chains',4,@(x)isnumeric(x)&&x>=1);
addParameter(p,'learn_nu',true,@islogical);
addParameter(p,'nu_step',0.25,@(x)isnumeric(x)&&x>0);
addParameter(p,'seed',6689,@isnumeric);

% 三件套
addParameter(p,'add_quadratic_idx',[],@(v)isnumeric(v)&&isvector(v));
addParameter(p,'add_spline_idx',[],@(v)isnumeric(v)&&isvector(v));
addParameter(p,'spline_knots',3,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'heteroskedastic_idx',[],@(v)isnumeric(v)&&isvector(v));
addParameter(p,'rand_slope_prov_idx',[],@(v)isnumeric(v)&&isvector(v));
addParameter(p,'rand_slope_var_idx',[],@(v)isnumeric(v)&&isvector(v));

% 先验强度 / 初始化
addParameter(p,'beta_sd',1.0,@(x)isnumeric(x)&&x>0);
addParameter(p,'alpha_sd',5.0,@(x)isnumeric(x)&&x>0);
addParameter(p,'ranef_scale',2.5,@(x)isnumeric(x)&&x>0);
addParameter(p,'init_nu',7,@(x)isnumeric(x)&&x>3);
addParameter(p,'init_jitter',0.2,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
parse(p,varargin{:});
args = p.Results;
rng(args.seed);

% ================== 载入缓存 ==================
cache_path = fullfile(args.results_root,'cache','hbm_cache.mat');
if ~isfile(cache_path), error('未找到缓存 %s（先运行 prep_hbm_inputs）', cache_path); end
Sdata = load(cache_path); cache = Sdata.cache;

out_root = fullfile(args.results_root,'M0_main');
if ~isfolder(out_root), mkdir(out_root); end

% ================== 基本量 ==================
N  = cache.N;
X0 = cache.X_fe;                 % FE（已标准化）
K0 = size(X0,2);
P  = cache.P; V = cache.V;
prov_id = cache.prov_id(:);
var_id  = cache.var_id(:);

% 连续协变量是 FE 前 Kc 列
Kc = numel(cache.cont_orig);
cont_cols = 1:Kc;

% ================== 扩展 FE（全样本） ==================
X = X0; fe_names = cache.active_fe_names;

% 二次项
qidx = intersect(unique(args.add_quadratic_idx(:)'), cont_cols);
for j = qidx
    z = X0(:,j).^2; mu = nanmean(z); sd = nanstd(z); if ~isfinite(sd)||sd==0, sd=1; end
    X = [X, (z-mu)/sd]; %#ok<AGROW>
    fe_names{end+1} = sprintf('%s^2', fe_names{j}); %#ok<AGROW>
end

% hinge 样条
sidx = intersect(unique(args.add_spline_idx(:)'), cont_cols);
KnotNum = args.spline_knots;
for j = sidx
    xj = X0(:,j);
    if KnotNum>0
        qs = linspace(0,1,KnotNum+2); qs = qs(2:end-1);
        knots = quantile(xj(isfinite(xj)), qs);
        for kk=1:numel(knots)
            z = max(0, xj - knots(kk));
            mu = nanmean(z); sd = nanstd(z); if ~isfinite(sd)||sd==0, sd=1; end
            X = [X, (z-mu)/sd]; %#ok<AGROW>
            fe_names{end+1} = sprintf('hinge(%s@q%.2g)', fe_names{j}, qs(kk)); %#ok<AGROW>
        end
    end
end
K = size(X,2);

% ================== 校验随机斜率/异方差索引 ==================
rs_prov_idx = intersect(unique(args.rand_slope_prov_idx(:)'), cont_cols);
rs_var_idx  = intersect(unique(args.rand_slope_var_idx(:)'),  cont_cols);
hs_idx_cont = intersect(unique(args.heteroskedastic_idx(:)'), cont_cols);

if ~isempty(rs_prov_idx), fprintf('[HBM] 随机斜率(省) 采用基础连续列索引: %s\n', mat2str(rs_prov_idx)); end
if ~isempty(hs_idx_cont), fprintf('[HBM] 异方差列 采用基础连续列索引: %s\n', mat2str(hs_idx_cont)); end

bad_rs = setdiff(unique(args.rand_slope_prov_idx(:)'), cont_cols);
bad_hs = setdiff(unique(args.heteroskedastic_idx(:)'),  cont_cols);
if ~isempty(bad_rs) || ~isempty(bad_hs)
    cont_list = cache.cont_orig(:);
    fprintf('>>> 连续基础列清单 (1..Kc):\n');
    for ii=1:numel(cont_list), fprintf('  %2d : %s\n', ii, cont_list{ii}); end
    error('rand_slope_prov_idx/heteroskedastic_idx 必须落在 1..Kc。错索引 rs=%s, hs=%s', ...
        mat2str(bad_rs), mat2str(bad_hs));
end
JP = numel(rs_prov_idx); JV = numel(rs_var_idx); H = numel(hs_idx_cont);

% 全样本索引（仅用于打印/校验）
% idx_prov_full = accumarray(prov_id, (1:N)', [P,1], @(v){v});
% idx_var_full  = accumarray(var_id,  (1:N)', [V,1], @(v){v});

targets = fieldnames(cache.Y);

% ================== 目标循环 ==================
for ti = 1:numel(targets)
    yname = targets{ti};
    y_tr  = cache.Y.(yname);

    ok = isfinite(y_tr) & all(isfinite(X),2);
    idx_ok = find(ok);
    N_use  = numel(idx_ok);
    if N_use < 20
        warning('[%s] 可用样本过少（N=%d），跳过。', yname, N_use);
        continue;
    end

    % 子样本作业面（同一坐标系）
    mu_y     = mean(y_tr(ok),'omitnan');     % 变换后尺度
    Xuse     = X(ok,:);                      % N_use × K
    X0use    = X0(ok,:);                     % N_use × K0（仅当你想把异方差基于基础列时备用）
    yuse     = y_tr(ok) - mu_y;              % N_use × 1
    prov_use = prov_id(ok);                  % N_use × 1
    var_use  = var_id(ok);                   % N_use × 1

    % 每个省/变异株的“行号”集合（严格在 1..N_use）
    idx_prov_use = cell(P,1);
    for p_i=1:P, idx_prov_use{p_i} = find(prov_use==p_i); end
    idx_var_use = cell(V,1);
    for v_i=1:V, idx_var_use{v_i} = find(var_use==v_i);  end

    % 异方差设计（全部使用 Xuse 坐标）
    Z_hs_use = [];
    if H>0
        Z_hs_use = Xuse(:, hs_idx_cont);
    end

    % 采样设置
    S_warm = args.num_warmup; S_keep = args.num_samples; C = args.chains;

    ydir = fullfile(out_root, yname);
    if ~isfolder(ydir), mkdir(ydir); end

    % ============= 多链采样 =============
    all_draws = cell(C,1);
    for ch = 1:C
        rng(args.seed + 1000*ti + ch);

        % 初始化
        alpha = 0 + args.init_jitter*randn();
        beta  = args.init_jitter*randn(K,1);
        b_prov = args.init_jitter*randn(P,1);
        b_var  = args.init_jitter*randn(V,1);
        sigma  = exp(0.1*randn());
        tau_prov = exp(0.1*randn());
        tau_var  = exp(0.1*randn());
        nu = args.init_nu;
        omega = ones(N_use,1);

        gamma_hs = zeros(H,1);
        b_prov_slope = zeros(P,JP); b_var_slope = zeros(V,JV);
        tau_prov_slope = ones(JP,1); tau_var_slope = ones(JV,1);

        % 先验
        prec_alpha = 1/args.alpha_sd^2;
        prec_beta  = ones(K,1) / args.beta_sd^2;
        Prec0      = diag([prec_alpha; prec_beta]);
        m0         = zeros(1+K,1);

        % 存储
        A_alpha = nan(S_keep,1); A_sigma = nan(S_keep,1);
        A_tau_p = nan(S_keep,1); A_tau_v = nan(S_keep,1);
        A_nu    = nan(S_keep,1);
        A_beta  = nan(S_keep, K);
        A_bprov = nan(S_keep, P); A_bvar = nan(S_keep, V);
        A_mu    = nan(S_keep, N);
        A_gamma = nan(S_keep, H);
        A_tau_ps = nan(S_keep, JP); A_tau_vs = nan(S_keep, JV);
        A_bprov_s = nan(S_keep, P*JP); A_bvar_s = nan(S_keep, V*JV);

        entered_keep = false;
        S_total = S_warm + S_keep;

        for it = 1:S_total
            if it==1
                fprintf('[%s][chain %d] S_warm=%d S_keep=%d | K=%d(base=%d, add=%d) H=%d JP=%d JV=%d\n', ...
                    yname, ch, S_warm, S_keep, K, K0, K-K0, H, JP, JV);
            elseif it==S_warm
                fprintf('[%s][chain %d] warmup end\n', yname, ch);
            elseif it==S_warm+1
                fprintf('[%s][chain %d] enter keep\n', yname, ch);
            end

            % ---- 0) 线性预测（全体）----
            mu_use = alpha + Xuse*beta + b_prov(prov_use) + b_var(var_use);
            if JP>0
                for j=1:JP
                    mu_use = mu_use + b_prov_slope(prov_use,j) .* Xuse(:, rs_prov_idx(j));
                end
            end
            if JV>0
                for j=1:JV
                    mu_use = mu_use + b_var_slope(var_use,j)  .* Xuse(:, rs_var_idx(j));
                end
            end
            resid2 = (yuse - mu_use).^2;

            % ---- 1) 异方差 ----
            if H>0
                zeta = Z_hs_use * gamma_hs;
                zeta(~isfinite(zeta)) = 0;
            else
                zeta = zeros(N_use,1);
            end

            % ---- 2) omega | ... ----
            a_post = (nu+1)/2;
            b_post = (nu + resid2 ./ max(sigma^2 .* exp(zeta), 1e-12))/2;
            omega  = gamrnd(a_post, 1./b_post);
            omega(~isfinite(omega) | omega<=0) = 1;

            % ---- 3) [alpha; beta] | ... ----
            W = omega(:) ./ max(sigma^2 .* exp(zeta), 1e-12);
            Xaug = [ones(N_use,1), Xuse];
            slope_full = slope_contrib_only(JP,JV,rs_prov_idx,rs_var_idx, ...
                            b_prov_slope,b_var_slope, Xuse, (1:N_use)', prov_use, var_use);
            y_til = yuse(:) - ( b_prov(prov_use) + b_var(var_use) + slope_full );
            WX = Xaug .* W;
            XtWX = Xaug' * WX;
            XtWy = Xaug' * (W .* y_til);
            Prec = XtWX + Prec0; rhs = XtWy + Prec0*m0;
            PP = (Prec + Prec')/2; jitter=1e-10; okchol=false;
            for ttry=1:5
                try
                    L = chol(PP + jitter*eye(1+K),'lower'); okchol=true; break;
                catch
                    jitter=jitter*10;
                end
            end
            if ~okchol, L = chol(eye(1+K),'lower'); rhs=zeros(1+K,1); end
            mu_ab = L'\(L\rhs); ab = mu_ab + L'\randn(1+K,1);
            alpha = ab(1); beta = ab(2:end);

            % ---- 4) b_prov | ... （子集安全）----
            for p_i=1:P
                idx = idx_prov_use{p_i};
                if isempty(idx)
                    b_prov(p_i) = randn();
                    continue;
                end
                % 子块
                X_loc   = Xuse(idx,:);  y_loc   = yuse(idx);
                pv_loc  = prov_use(idx); vv_loc = var_use(idx);
                W_loc   = W(idx);
                slope_l = slope_contrib_only(JP,JV,rs_prov_idx,rs_var_idx, ...
                            b_prov_slope,b_var_slope, X_loc, pv_loc, vv_loc); % 旧签名：Xsub,pid,vid
                mu0 = y_loc - alpha - X_loc*beta - b_var(vv_loc) - slope_l;
                varp  = 1/(1/max(tau_prov^2,1e-12) + sum(W_loc));
                meanp = varp * sum(W_loc .* mu0);
                b_prov(p_i) = meanp + sqrt(varp)*randn();
            end

            % ---- 5) b_var | ... （子集安全）----
            for v_i=1:V
                idx = idx_var_use{v_i};
                if isempty(idx)
                    b_var(v_i) = randn();
                    continue;
                end
                X_loc   = Xuse(idx,:);  y_loc   = yuse(idx);
                pv_loc  = prov_use(idx); vv_loc = var_use(idx);
                W_loc   = W(idx);
                slope_l = slope_contrib_only(JP,JV,rs_prov_idx,rs_var_idx, ...
                            b_prov_slope,b_var_slope, X_loc, pv_loc, vv_loc);
                mu0 = y_loc - alpha - X_loc*beta - b_prov(pv_loc) - slope_l;
                varv  = 1/(1/max(tau_var^2,1e-12) + sum(W_loc));
                meanv = varv * sum(W_loc .* mu0);
                b_var(v_i) = meanv + sqrt(varv)*randn();
            end

            % ---- 6) 随机斜率（省）（子集安全）----
            for j=1:JP
                col = rs_prov_idx(j);
                for p_i=1:P
                    idx = idx_prov_use{p_i};
                    if isempty(idx)
                        b_prov_slope(p_i,j) = randn()*tau_prov_slope(j);
                        continue;
                    end
                    X_loc   = Xuse(idx,:);  y_loc   = yuse(idx);
                    pv_loc  = prov_use(idx); vv_loc = var_use(idx);
                    W_loc   = W(idx);       xj_loc = X_loc(:, col);
                    % slope（扣除“当前 j”的项后）——用“先加后减”的办法避免写两套函数
                    slope_l = slope_contrib_only(JP,JV,rs_prov_idx,rs_var_idx, ...
                                 b_prov_slope,b_var_slope, X_loc, pv_loc, vv_loc);
                    slope_l = slope_l - b_prov_slope(p_i,j) * xj_loc .* (pv_loc==p_i);
                    mu0 = y_loc - alpha - X_loc*beta - b_prov(pv_loc) - b_var(vv_loc) - slope_l;
                    prec = 1/max(tau_prov_slope(j)^2,1e-12) + sum(W_loc .* (xj_loc.^2));
                    varp = 1/max(prec,1e-12);
                    meanp= varp * sum(W_loc .* xj_loc .* mu0);
                    b_prov_slope(p_i,j) = meanp + sqrt(varp)*randn();
                end
            end

            % ---- 7) 随机斜率（变异株）（子集安全）----
            for j=1:JV
                col = rs_var_idx(j);
                for v_i=1:V
                    idx = idx_var_use{v_i};
                    if isempty(idx)
                        b_var_slope(v_i,j) = randn()*tau_var_slope(j);
                        continue;
                    end
                    X_loc   = Xuse(idx,:);  y_loc   = yuse(idx);
                    pv_loc  = prov_use(idx); vv_loc = var_use(idx);
                    W_loc   = W(idx);       xj_loc = X_loc(:, col);
                    slope_l = slope_contrib_only(JP,JV,rs_prov_idx,rs_var_idx, ...
                                 b_prov_slope,b_var_slope, X_loc, pv_loc, vv_loc);
                    slope_l = slope_l - b_var_slope(v_i,j) * xj_loc .* (vv_loc==v_i);
                    mu0 = y_loc - alpha - X_loc*beta - b_prov(pv_loc) - b_var(vv_loc) - slope_l;
                    prec = 1/max(tau_var_slope(j)^2,1e-12) + sum(W_loc .* (xj_loc.^2));
                    varp = 1/max(prec,1e-12);
                    meanv= varp * sum(W_loc .* xj_loc .* mu0);
                    b_var_slope(v_i,j) = meanv + sqrt(varv)*randn();
                end
            end

            % ---- 8) sigma（异方差）----
            sigma = slice_sample_pos(@(s) logpost_sigma_het(s, resid2, omega, zeta, args.ranef_scale), sigma);
            if ~(isfinite(sigma)&&sigma>0), sigma=1.0; end

            % ---- 9) tau_prov / tau_var ----
            tau_prov = slice_sample_pos(@(t) logpost_tau(t, b_prov, args.ranef_scale), tau_prov);
            if ~(isfinite(tau_prov)&&tau_prov>0), tau_prov=1.0; end
            tau_var  = slice_sample_pos(@(t) logpost_tau(t, b_var,  args.ranef_scale), tau_var);
            if ~(isfinite(tau_var)&&tau_var>0),  tau_var=1.0;  end

            % ---- 10) 随机斜率尺度 ----
            for j2=1:JP
                tau_prov_slope(j2) = slice_sample_pos(@(t) logpost_tau(t, b_prov_slope(:,j2), args.ranef_scale), tau_prov_slope(j2));
                if ~(isfinite(tau_prov_slope(j2))&&tau_prov_slope(j2)>0), tau_prov_slope(j2)=1.0; end
            end
            for j2=1:JV
                tau_var_slope(j2)  = slice_sample_pos(@(t) logpost_tau(t, b_var_slope(:,j2),  args.ranef_scale), tau_var_slope(j2));
                if ~(isfinite(tau_var_slope(j2))&&tau_var_slope(j2)>0),  tau_var_slope(j2)=1.0;  end
            end

            % ---- 11) 异方差系数 gamma（全实轴）----
            if H>0
                for j2=1:H
                    gamma_hs(j2) = slice_sample_real(@(g) logpost_gamma_j(g, j2, gamma_hs, Z_hs_use, resid2, omega, sigma), gamma_hs(j2));
                end
            end

            % ---- 12) nu（可选）----
            if args.learn_nu
                nu = mh_update_nu(nu, omega, args.nu_step);
                if ~(isfinite(nu)&&nu>3), nu=7; end
            end

            % ---- 保存（仅保留期）----
            if it > S_warm
                k = it - S_warm; entered_keep = true;
                A_alpha(k)=alpha; A_sigma(k)=sigma; A_tau_p(k)=tau_prov; A_tau_v(k)=tau_var; A_nu(k)=nu;
                A_beta(k,:)=beta(:).'; A_bprov(k,:)=b_prov(:).'; A_bvar(k,:)=b_var(:).';
                if H>0, A_gamma(k,:)=gamma_hs(:).'; end
                if JP>0, A_tau_ps(k,:)=tau_prov_slope(:).'; A_bprov_s(k,:)=reshape(b_prov_slope.',1,[]); end
                if JV>0, A_tau_vs(k,:)=tau_var_slope(:).';  A_bvar_s(k,:)=reshape(b_var_slope.',1,[]);  end

                mu_full = nan(N,1);
                mu_tmp = alpha + Xuse*beta + b_prov(prov_use) + b_var(var_use);
                if JP>0
                    for j2=1:JP
                        mu_tmp = mu_tmp + b_prov_slope(prov_use,j2) .* Xuse(:, rs_prov_idx(j2));
                    end
                end
                if JV>0
                    for j2=1:JV
                        mu_tmp = mu_tmp + b_var_slope(var_use,j2)  .* Xuse(:, rs_var_idx(j2));
                    end
                end
                mu_tmp = mu_tmp + mu_y;
                mu_full(idx_ok) = mu_tmp;
                A_mu(k,:) = mu_full(:).';
                if k==1
                    fprintf('[%s][chain %d] kept#1: a=%.3f, s=%.3f, tp=%.3f, tv=%.3f, nu=%.2f | H=%d JP=%d JV=%d\n', ...
                        yname, ch, A_alpha(k),A_sigma(k),A_tau_p(k),A_tau_v(k),A_nu(k), H,JP,JV);
                end
            end
        end % it

        if ~entered_keep, warning('[%s][chain %d] 未进入保留阶段。', yname, ch); end

        % 输出 table
        vars = struct();
        vars.alpha=A_alpha; vars.sigma=A_sigma; vars.sigma_prov=A_tau_p; vars.sigma_var=A_tau_v; vars.nu=A_nu;
        for j=1:K, vars.(sprintf('beta_%d',j))=A_beta(:,j); end
        for j=1:P, vars.(sprintf('b_prov_%d',j))=A_bprov(:,j); end
        for j=1:V, vars.(sprintf('b_var_%d',j)) =A_bvar(:,j); end
        for n=1:N, vars.(sprintf('mu_%d',n))    =A_mu(:,n);   end
        if H>0, for j=1:H,  vars.(sprintf('gamma_hs_%d',j)) = A_gamma(:,j); end, end
        if JP>0
            for j=1:JP, vars.(sprintf('tau_prov_slope_%d',j))=A_tau_ps(:,j); end
            for j=1:JP, for p_i=1:P
                col=(j-1)*P+p_i; vars.(sprintf('b_prov_slope_%d_%d',j,p_i))=A_bprov_s(:,col);
            end, end
        end
        if JV>0
            for j=1:JV, vars.(sprintf('tau_var_slope_%d',j))=A_tau_vs(:,j); end
            for j=1:JV, for v_i=1:V
                col=(j-1)*V+v_i; vars.(sprintf('b_var_slope_%d_%d',j,v_i))=A_bvar_s(:,col);
            end, end
        end
        vars.chain = ch * ones(S_keep,1);
        vars.iter  = (1:S_keep)';

        Tch = table(); fns = fieldnames(vars);
        for f=1:numel(fns), Tch.(fns{f}) = vars.(fns{f}); end
        all_draws{ch} = Tch;
        fprintf('[MAT-HBM+] %s chain %d done. kept=%d\n', yname, ch, S_keep);
    end % chain

    % 合并多链并保存
    draws = vertcat(all_draws{:});
    save(fullfile(ydir,'draws_raw.mat'),'draws','-v7.3');

    % 记录扩展特征名
    meta = struct();
    meta.fe_names = fe_names; meta.K=K; meta.K_base=K0; meta.Kc=Kc;
    meta.quadratic_idx = qidx; meta.spline_idx = sidx;
    meta.heteroskedastic_idx = hs_idx_cont;
    meta.rand_slope_prov_idx = rs_prov_idx;
    meta.rand_slope_var_idx  = rs_var_idx;
    save(fullfile(ydir,'fe_meta.mat'),'meta');

    fid = fopen(fullfile(ydir,'fit_done.txt'),'w'); fprintf(fid,'ok\n'); fclose(fid);
end

fprintf('[fit_hbm_main_m0_matlab+] 全部完成。\n');

% ================== 辅助函数 ==================

function v = slope_contrib_only(JP,JV,rsP,rsV,bpS,bvS,Xblk,varargin)
% 两种用法：
%   (1) slope_contrib_only(..., Xuse, row_idx, pid, vid)
%   (2) slope_contrib_only(..., Xsub,        pid, vid)
    if numel(varargin) == 3
        row_idx = varargin{1}; pid = varargin{2}; vid = varargin{3};
        Xsub = Xblk(row_idx, :);
    elseif numel(varargin) == 2
        pid = varargin{1}; vid = varargin{2};
        Xsub = Xblk;
    else
        error('slope_contrib_only:arity','参数需为 (Xuse,row_idx,pid,vid) 或 (Xsub,pid,vid)。');
    end
    pid = pid(:); vid = vid(:);
    n = size(Xsub,1);
    if numel(pid) ~= n || numel(vid) ~= n
        error('slope_contrib_only:lenMismatch','pid/vid 与 Xsub 行数不一致');
    end
    v = zeros(n,1);
    if JP > 0
        for j = 1:JP
            col = rsP(j);
            if col<1 || col>size(Xsub,2), error('省随机斜率列越界'); end
            v = v + bpS(pid,j) .* Xsub(:,col);
        end
    end
    if JV > 0
        for j = 1:JV
            col = rsV(j);
            if col<1 || col>size(Xsub,2), error('株随机斜率列越界'); end
            v = v + bvS(vid,j) .* Xsub(:,col);
        end
    end
end

function lp = logpost_sigma_het(s, resid2, omega, zeta, scale_h)
    if ~(isfinite(s) && s>0), lp = -inf; return; end
    s2 = max(s^2,1e-12);
    ll = -0.5*sum( log(s2) + zeta - log(max(omega,1e-300)) + (omega .* resid2) ./ (s2 .* exp(zeta)) );
    nu0 = 3; sc = scale_h;
    lp_prior = - (nu0+1)/2 * log(1 + (s/sc)^2/nu0);
    lp = ll + lp_prior;
end

function lp = logpost_tau(t, bvec, scale_h)
    if ~(isfinite(t) && t>0), lp = -inf; return; end
    ll = -numel(bvec)*log(max(t,1e-12)) - sum(bvec.^2) / (2*max(t^2,1e-12));
    nu0 = 3; sc = scale_h;
    lp_prior = - (nu0+1)/2 * log(1 + (t/sc)^2/nu0);
    lp = ll + lp_prior;
end

function lp = logpost_gamma_j(g, j, gamma_cur, Z, resid2, omega, sigma)
    if ~isfinite(g), lp = -inf; return; end
    gamma_tmp = gamma_cur; gamma_tmp(j) = g;
    zeta = Z * gamma_tmp;
    s2 = max(sigma^2,1e-12);
    ll = -0.5*sum( zeta - log(max(omega,1e-300)) + (omega .* resid2) ./ (s2 .* exp(zeta)) );
    lp_prior = -0.5 * sum(gamma_tmp.^2);
    lp = ll + lp_prior;
end

function x_new = slice_sample_pos(logf, x0)
    if ~(isfinite(x0) && x0>0), x0 = 1; end
    w=1; m=20; maxit=100;
    ly = logf(x0); if ~isfinite(ly), ly = -1e6; end
    log_y = ly + log(rand()); u=rand();
    L = max(x0 - w*u, eps); R = x0 + w*(1-u);
    J = floor(m*rand()); Kk = (m-1)-J;
    while J>0 && logf(L) > log_y, L = max(L - w, eps); J = J - 1; end
    while Kk>0 && logf(R) > log_y, R = R + w;           Kk = Kk - 1; end
    for ssit=1:maxit
        x1 = L + rand()*(R-L); lf = logf(x1);
        if ~isfinite(lf), if x1>x0, R=x1; else, L=x1; end; continue; end
        if lf >= log_y, x_new = x1; if ~(isfinite(x_new)&&x_new>0), x_new=x0; end; return;
        else, if x1>x0, R=x1; else, L=x1; end
        end
    end
    x_new = max(x0,1e-6);
end

function nu_new = mh_update_nu(nu_old, omega, step)
    if ~(isfinite(nu_old)&&nu_old>3), nu_old=7; end
    z_old = log(nu_old - 3); z_new = z_old + step*randn();
    nu_cand = 3 + exp(z_new);
    if ~(isfinite(nu_cand)&&nu_cand>3), nu_new=nu_old; return; end
    logpost = @(nu) logp_omega_given_nu(omega, nu) + logprior_nu(nu);
    a = (logpost(nu_cand) + log(nu_cand-3)) - (logpost(nu_old) + log(nu_old-3));
    if log(rand()) < a, nu_new = nu_cand; else, nu_new = nu_old; end
end
function lp = logp_omega_given_nu(omega, nu)
    a = nu/2; b = nu/2;
    lp = numel(omega)*( a*log(b) - gammaln(a) ) + sum( (a-1)*log(max(omega,1e-300)) - b*omega );
end
function lp = logprior_nu(nu)
    x = nu - 3; if x<=0 || ~isfinite(x), lp=-inf; return; end
    k=2; r=0.1; lp = (k-1)*log(x) - r*x - gammaln(k) + k*log(r);
end

function x_new = slice_sample_real(logf, x0)
    if ~isfinite(x0), x0 = 0; end
    w=0.5; m=10; maxit=50;
    fx0 = logf(x0); if ~isfinite(fx0), fx0 = -1e-6; end
    log_y = fx0 + log(rand()); u=rand();
    L = x0 - w*u; R = x0 + w*(1-u);
    J = floor(m*rand()); Kk = (m-1)-J;
    while J>0 && logf(L) > log_y, L = L - w; J = J - 1; end
    while Kk>0 && logf(R) > log_y, R = R + w; Kk = Kk - 1; end
    for ssit=1:maxit
        x1 = L + rand()*(R-L); lf = logf(x1);
        if isfinite(lf) && lf >= log_y, x_new = x1; return; end
        if x1 > x0, R = x1; else, L = x1; end
    end
    x_new = x0;
end
end
