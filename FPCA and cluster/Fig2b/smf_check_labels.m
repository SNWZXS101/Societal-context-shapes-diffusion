function smf_check_labels(matfile, Kwant)
S = load(matfile);
flds = fieldnames(S);
if numel(flds)==1 && isstruct(S.(flds{1})), S = S.(flds{1}); end
get = @(T,names) local_pick(T, string(names));

F = get(S,["F_norm","F","curves","X","Y"]);
u = get(S,["u","grid","time_grid","U"]);
if ~isempty(F)&&~isempty(u) && size(F,1)~=numel(u) && size(F,2)==numel(u), F=F.'; end
N = size(F,2);

clu = []; if isstruct(S) && isfield(S,'cluster'), clu=S.cluster; end
labC = []; KC=[];
if ~isempty(clu)
    labC = get(clu,["labels","idx","cluster_id","assignments"]); 
    if isfield(clu,'K'), KC=clu.K; else, KC = max(labC); end
end

labB = []; KB=[]; 
if isfield(S,'byK') && ~isempty(S.byK)
    kk = arrayfun(@(z) z.K, S.byK);
    p = find(kk==Kwant,1);
    if ~isempty(p), labB = S.byK(p).labels; KB = Kwant; end
end

fprintf('F_norm N = %d\n', N);
if ~isempty(labC)
    fprintf('cluster.labels N = %d (K=%d) | counts: %s\n', numel(labC), KC, mat2str(histcounts(labC,0.5:1:KC+0.5)));
else
    fprintf('cluster.labels: not found\n');
end
if ~isempty(labB)
    fprintf('byK(K=%d).labels N = %d | counts: %s\n', KB, numel(labB), mat2str(histcounts(labB,0.5:1:KB+0.5)));
else
    fprintf('byK(K=%d): not found\n', Kwant);
end

if ~isempty(labC) && ~isempty(labB)
    if numel(labC)~=numel(labB)
        warning('cluster 和 byK(K=%d) 的样本总数不同！这就是你看到分布变动的原因。', Kwant);
    end
    if numel(labC)~=N
        warning('cluster.labels (N=%d) 与 F_norm 列数 (N=%d) 不同。', numel(labC), N);
    end
    if numel(labB)~=N
        warning('byK(K=%d).labels (N=%d) 与 F_norm 列数 (N=%d) 不同。', Kwant, numel(labB), N);
    end
end
end

function val = local_pick(T, keys)
val=[]; if isempty(T), return; end
if istable(T)
    for i=1:numel(keys), k=char(keys(i)); if any(strcmpi(T.Properties.VariableNames,k)), val=T.(k); return; end, end
elseif isstruct(T)
    f=fieldnames(T);
    for i=1:numel(keys), k=char(keys(i)); hit=f(strcmpi(f,k)); if ~isempty(hit), val=T.(hit{1}); return; end, end
end
end
