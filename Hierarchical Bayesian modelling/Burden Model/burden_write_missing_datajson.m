function ok = burden_write_missing_datajson(results_root, target_dir)
% 目的：为 results_root/burden/target_dir/ 补写缺失的 data.json
% 仅依赖：
%   - results_root/cache/hbm_cache.mat  （提供 prov_id / var_id / X_fe）
%   - results_root/burden/target_dir/ppc.csv （提供 N，用于对齐样本）
%
% 写出内容：
%   { "N": <int>, "province_id": [..], "variant_id": [..] }

out_dir  = fullfile(results_root,'burden',target_dir);
json_path = fullfile(out_dir,'data.json');
if isfile(json_path), ok=true; return; end

% 1) 读 cache
cache_path = fullfile(results_root,'cache','hbm_cache.mat');
if ~isfile(cache_path)
    warning('[repair] 缺少 %s，无法重建 data.json。', cache_path); ok=false; return;
end
C = load(cache_path); cache = C.cache;

% 2) 读 ppc.csv（拿到 N）
ppc_path = fullfile(out_dir,'ppc.csv');
if ~isfile(ppc_path)
    warning('[repair] 缺少 %s，无法推断 N。', ppc_path); ok=false; return;
end
PPC = readtable(ppc_path);
Nppc = height(PPC);

% 3) 用 X_fe 的有限性作为拟合有效样本掩码（与拟合时一致）
ok_mask = all(isfinite(cache.X_fe), 2);
idx_all = find(ok_mask);

if numel(idx_all) < Nppc
    warning('[repair] 可用行数(%d) < ppc 行数(%d)，放弃重建。', numel(idx_all), Nppc);
    ok=false; return;
end

% 与拟合顺序一致地取前 Nppc 行
idx_use = idx_all(1:Nppc);

data_struct = struct();
data_struct.N = int32(numel(idx_use));
data_struct.province_id = int32(cache.prov_id(idx_use));
data_struct.variant_id  = int32(cache.var_id(idx_use));

try
    fid = fopen(json_path,'w');
    fwrite(fid, jsonencode(data_struct,'PrettyPrint',true));
    fclose(fid);
    ok = true;
    fprintf('[repair] 已重建 %s（N=%d）\n', json_path, data_struct.N);
catch
    warning('[repair] 写 %s 失败。', json_path); ok=false;
end
end
