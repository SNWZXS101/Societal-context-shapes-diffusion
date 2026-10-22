function key = eng_name_to_key(eng_name)
s = lower(string(eng_name));
% 特例先判
if contains(s,'inner mongolia'), key="InnerMongolia"; return; end
if contains(s,'tibet'),          key="Tibet";         return; end
if contains(s,'hong kong'),      key="HongKong";      return; end
if contains(s,'macau') || contains(s,'macao')
    key="Macau"; return;
end
% 去常见后缀
rm = {' province',' autonomous region',' municipality',' special administrative region', ...
      ' uygur',' hui',' zhuang'};
for i=1:numel(rm), s = replace(s, rm{i}, ''); end
s = regexprep(s,'\s+','');
switch s
    case 'beijing',key="Beijing"; case 'tianjin',key="Tianjin"; case 'hebei',key="Hebei";
    case 'shanxi',key="Shanxi";   case 'shaanxi',key="Shaanxi";
    case 'liaoning',key="Liaoning"; case 'jilin',key="Jilin"; case 'heilongjiang',key="Heilongjiang";
    case 'shanghai',key="Shanghai"; case 'jiangsu',key="Jiangsu"; case 'zhejiang',key="Zhejiang";
    case 'anhui',key="Anhui"; case 'fujian',key="Fujian"; case 'jiangxi',key="Jiangxi";
    case 'shandong',key="Shandong"; case 'henan',key="Henan"; case 'hubei',key="Hubei";
    case 'hunan',key="Hunan"; case 'guangdong',key="Guangdong"; case 'guangxi',key="Guangxi";
    case 'hainan',key="Hainan"; case 'chongqing',key="Chongqing"; case 'sichuan',key="Sichuan";
    case 'guizhou',key="Guizhou"; case 'yunnan',key="Yunnan"; case 'gansu',key="Gansu";
    case 'qinghai',key="Qinghai"; case 'ningxia',key="Ningxia"; case 'xinjiang',key="Xinjiang";
    case 'taiwan',key="Taiwan";
    otherwise, key="";
end
end
