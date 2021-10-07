%% Export values from porated and nonporated locations into tables to be used for Machine Learning
clear all; clc
mkdir DataTables

varnameA = {'parea','thick','mcurv','dipln','charg','P2_cosTAmean'};
varnameB = {'FS','MU','PU','CHOL','PC','PE','SM','GM','CE','LPC','DAG','PS','PI','PA','PIP'};

PM_list = {'APM-dep'}; sim = 'equil10nsBeforeEField'; rr = 15; frames0 = 0:20;

for ii = 1:length(PM_list)
    PM = PM_list{ii};
    load([PM,'/forML_',sim,'_rr',num2str(rr),'.mat'])
    load(['PoreLocations/',PM,'.mat'])
    
    for mem = 1%:4
        xp = poreloci_rel{mem}(:,1);
        yp = poreloci_rel{mem}(:,2);
        tp = poreloci{mem}(:,2);
        xnp = noporeloci_rel{mem}(:,1);
        ynp = noporeloci_rel{mem}(:,2);
        tnp = zeros(size(xnp)) + nan;
        
        Xp = repmat(xp,length(frames0),1);
        Yp = repmat(yp,length(frames0),1);
        Tp = repmat(tp,length(frames0),1);
        Fp = sort(repmat((1:length(frames0))',length(xp),1));
        Xnp = repmat(xnp,length(frames0),1);
        Ynp = repmat(ynp,length(frames0),1);
        Tnp = repmat(tnp,length(frames0),1);
        Fnp = sort(repmat((1:length(frames0))',length(xnp),1));
       
        data_porated = [Xp, Yp, Tp, Fp];
        data_nonporated = [Xnp, Ynp, Tnp, Fnp];
        header = {'x_rel','y_rel','tpore','frame'};

        % General properties: Mean of leaflets
        varname = varnameA;
        for j = 1:length(varname)
            data_porated = [data_porated, poratedA{mem}.mean.(varname{j})];
            data_nonporated = [data_nonporated, nonporatedA{mem}.mean.(varname{j})];
            header = [header, [varname{j},'_mean']];
        end
        % Lipid densities
        varname = varnameB;
        for j = 1:length(varname)
            data_porated = [data_porated, poratedB{mem}.sum.(varname{j})];
            data_nonporated = [data_nonporated, nonporatedB{mem}.sum.(varname{j})];
            header = [header, varname{j}];
        end
                      
        data_porated_table = array2table(data_porated,'VariableNames',header);
        data_nonporated_table = array2table(data_nonporated,'VariableNames',header);
        writetable(data_porated_table,['DataTables/',PM,'_mem',num2str(mem),'_porated','.csv'])
        writetable(data_nonporated_table,['DataTables/',PM,'_mem',num2str(mem),'_nonporated','.csv'])        
    end
end
