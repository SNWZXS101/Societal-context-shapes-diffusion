function check_feature_coverage(results_root)
if nargin<1, results_root='results_hbm'; end
S = load(fullfile(results_root,'cache','hbm_cache.mat')); C = S.cache;

fprintf('K=%d fixed effects:\n', size(C.X_fe,2));
for i=1:numel(C.fe_names_orig)
    fprintf('  #%02d %s\n', i, C.fe_names_orig{i});
end

badrow = ~all(isfinite(C.X_fe),2);
fprintf('行中含非有限FE的条目数：%d / %d\n', sum(badrow), C.N);

% 哪些列还有非有限值
badcols = ~all(isfinite(C.X_fe),1);
if any(badcols)
    fprintf('仍含非有限的FE列（列名: 非有限计数）:\n');
    for j=find(badcols)
        cnt = sum(~isfinite(C.X_fe(:,j)));
        fprintf('  - %s: %d\n', C.fe_names_orig{j}, cnt);
    end
else
    fprintf('所有FE列均为有限数。\n');
end

targets = fieldnames(C.Y);
for i=1:numel(targets)
    yname = targets{i};
    y = C.Y.(yname);
    y_nonmiss = sum(isfinite(y));
    ok = isfinite(y) & all(isfinite(C.X_fe),2);
    fprintf('[%s] y非NaN=%d / %d,  N_use=%d / %d\n', yname, y_nonmiss, C.N, sum(ok), C.N);
end
end
