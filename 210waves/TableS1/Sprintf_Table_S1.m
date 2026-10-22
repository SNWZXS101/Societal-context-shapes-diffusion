% 默认剔除三大波，保留 1 位小数
summarize_table_S1('Table_S1.xlsx');

% 不剔除三大波，保留 0 位小数
summarize_table_S1('Table_S1.xlsx','ExcludeMega',false,'Decimals',1);

% 只打印句子（静默数量信息）
summarize_table_S1('Table_S1.xlsx','Quiet',true);
