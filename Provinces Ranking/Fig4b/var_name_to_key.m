function key = var_name_to_key(var_name, eng_name)
% VAR_NAME（带声调拼音）→ 项目拼音键
% - 去声调（āáǎà→a 等）、去空格/连字符/点
% - 去后缀：sheng/shi/zizhiqu/tebiexingzhengqu/zhuangzu/huizu/weiwuer
% - 处理 Aomen/Macao、Xianggang/HongKong、Neimenggu/InnerMongolia、Xizang/Tibet
% - 对“shanxi”用 ENG_NAME 判别山西/陕西

s = string(var_name);
e = string(eng_name);

% 1) 统一小写、去空白与标点
s = lower(s);
s = regexprep(s, '[\s\-\._]+', '');

% 2) 去声调
rep = {
    'ā','a'; 'á','a'; 'ǎ','a'; 'à','a'; 'ā','a'; 'á','a'; 'ǎ','a'; 'à','a'; 'ɑ','a';
    'ē','e'; 'é','e'; 'ě','e'; 'è','e'; 'ē','e'; 'é','e'; 'ě','e'; 'è','e';
    'ī','i'; 'í','i'; 'ǐ','i'; 'ì','i'; 'ī','i'; 'í','i'; 'ǐ','i'; 'ì','i';
    'ō','o'; 'ó','o'; 'ǒ','o'; 'ò','o'; 'ō','o'; 'ó','o'; 'ǒ','o'; 'ò','o';
    'ū','u'; 'ú','u'; 'ǔ','u'; 'ù','u'; 'ū','u'; 'ú','u'; 'ǔ','u'; 'ù','u';
    'ǖ','u'; 'ǘ','u'; 'ǚ','u'; 'ǜ','u'; 'ü','u';
    'ń','n'; 'ň','n'; 'ǹ','n'
};
for i=1:size(rep,1), s = replace(s, rep{i,1}, rep{i,2}); end

% 3) 去后缀/族名
rm = {'sheng','shi','zizhiqu','tebiexingzhengqu','zhuangzu','huizu','weiwuer','uygur'};
for i=1:numel(rm), s = replace(s, rm{i}, ''); end

% 4) 特例归一
s = replace(s, "aomen", "macau");
s = replace(s, "macao", "macau");
s = replace(s, "xianggang", "hongkong");
s = replace(s, "neimenggu", "innermongolia");
s = replace(s, "xizang", "tibet");
s = replace(s, "ningxiahui", "ningxia");
s = replace(s, "ningxiahui zu", "ningxia"); % 保险
s = replace(s, "xinjiangweiwuer", "xinjiang");

% 5) 映射
switch s
    case 'beijing',        key = "Beijing";
    case 'tianjin',        key = "Tianjin";
    case 'hebei',          key = "Hebei";
    case 'shanxi'
        % 用英文名 disambiguate：Shaanxi vs Shanxi
        e2 = lower(e);
        if contains(e2,'shaanxi'), key="Shaanxi"; else, key="Shanxi"; end
    case 'innermongolia',  key = "InnerMongolia";
    case 'liaoning',       key = "Liaoning";
    case 'jilin',          key = "Jilin";
    case 'heilongjiang',   key = "Heilongjiang";
    case 'shanghai',       key = "Shanghai";
    case 'jiangsu',        key = "Jiangsu";
    case 'zhejiang',       key = "Zhejiang";
    case 'anhui',          key = "Anhui";
    case 'fujian',         key = "Fujian";
    case 'jiangxi',        key = "Jiangxi";
    case 'shandong',       key = "Shandong";
    case 'henan',          key = "Henan";
    case 'hubei',          key = "Hubei";
    case 'hunan',          key = "Hunan";
    case 'guangdong',      key = "Guangdong";
    case 'guangxi',        key = "Guangxi";
    case 'hainan',         key = "Hainan";
    case 'chongqing',      key = "Chongqing";
    case 'sichuan',        key = "Sichuan";
    case 'guizhou',        key = "Guizhou";
    case 'yunnan',         key = "Yunnan";
    case 'tibet',          key = "Tibet";
    case 'gansu',          key = "Gansu";
    case 'qinghai',        key = "Qinghai";
    case 'ningxia',        key = "Ningxia";
    case 'xinjiang',       key = "Xinjiang";
    case 'taiwan',         key = "Taiwan";
    case 'hongkong',       key = "HongKong";
    case 'macau',          key = "Macau";
    otherwise,             key = "";
end
end
