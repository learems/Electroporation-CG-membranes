function [porated_out,nonporated_out] = groupData(porated,nonporated,varname,memids)

tpbt_list = {'bt','tp','sum','mean','diff'};

for k = 1:length(tpbt_list)
    tpbt = tpbt_list{k};
    for j = 1:length(varname)
        porated_out.(tpbt).(varname{j}) = [];
        nonporated_out.(tpbt).(varname{j}) = [];
        for mem = memids
            porated_out.(tpbt).(varname{j}) = [porated_out.(tpbt).(varname{j}); porated{mem}.(tpbt).(varname{j})];
            nonporated_out.(tpbt).(varname{j}) = [nonporated_out.(tpbt).(varname{j}); nonporated{mem}.(tpbt).(varname{j})];
        end
    end
end
    

