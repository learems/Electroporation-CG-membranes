%% Extract local properties from porated and nonporated locations
addpath('./MatlabFunctions')

%% General
PM_list = {'APM-dep'}; sim = 'equil10nsBeforeEField'; frames = 0:20; 
rr = 15; % sigma of the gaussian kernel for smoothing

%% Analysis of equil10nsBeforeEField
varnameA = {'parea','thick','mcurv','dipln','charg','P2i_mean','P2_cosTAmean'};
varnameB = {'PC','PE','SM','GM','GM1','GM3','GMS','CE','LPC','CHOL','DAG','PS','PI','PA','PIP','FS','MU','PU1','PU2','PU'};

for ii = 1:length(PM_list)
    PM = PM_list{ii};
    
    % Load data on pore loctions
    load(['./PoreLocations/',PM,'.mat'])
    
    % Extract the values of selected variables from porated and nonporated locations
    clear porated1 nonporated1 porated2 nonporated2
    for mem = 1%:4
        load([PM,'/mem',num2str(mem),'/',sim,'/memsurf.mat'])
        % Values of variables A (varnameA) in porated and nonporated locations
        [poratedA{mem}, nonporatedA{mem}] = extractLocalValues(boxsize,memsurfb,memsurft,varnameA,frames,poreloci_rel{mem},noporeloci_rel{mem},'gaussKernel',rr);
        % Values of variables B (varnameB) in porated and nonporated locations
        [poratedB{mem}, nonporatedB{mem}] = extractLocalValues(boxsize,memsurfb,memsurft,varnameB,frames,poreloci_rel{mem},noporeloci_rel{mem},'gaussKernel',rr);
    end
    % Group data from all four membranes 
    % NOTE: This example considers a single membrane
%     [poratedA{5}, nonporatedA{5}] = groupData(poratedA,nonporatedA,varnameA,1:4);
%     [poratedB{5}, nonporatedB{5}] = groupData(poratedB,nonporatedB,varnameB,1:4);
    [poratedA{5}, nonporatedA{5}] = groupData(poratedA,nonporatedA,varnameA,1);
    [poratedB{5}, nonporatedB{5}] = groupData(poratedB,nonporatedB,varnameB,1);
    % Save
    save([PM,'/forML_',sim,'_rr',num2str(rr),'.mat'],'poratedA','nonporatedA','poratedB','nonporatedB','rr','varnameA','varnameB')    
end
