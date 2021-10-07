%% Load data and save data for analysis
% This script loads the data from txt files (extracted with MemSurfer) into mat files. 

%% Define what to analyse
clear all
PM_list = {'APM-dep'}; sim_list = {'equil10nsBeforeEField'};   frames = 0:20;

%%
% Variables extracted by MemSurfer: 
varname = {'verti','pnorm','parea','mcurv','thick','order','dipln','charg','resnum'};
% verti = coordinates of the vertices; pnorm = local vector of the membrane normal, local parea = area per lipid, mcurv = local mean curvature; 
% thick = local thickness; dipln = local cosine of the dipole angle; charg = local no. elementary charges; 
% resnum = residue number at a given vertex 

% List of all lipids
lipids = {'POPX', 'PIPX', 'DPPX', 'DAPC', 'DOPC', 'DPPC', 'OIPC', 'OUPC', 'PAPC', 'PEPC', 'PFPC', 'PIPC', 'POPC', 'PUPC', ...
          'DAPE', 'DOPE', 'DUPE', 'OAPE', 'OIPE', 'OUPE', 'PAPE', 'PIPE', 'POPE', 'PQPE', 'PUPE', ...
          'BNSM', 'DBSM', 'DPSM', 'DXSM', 'PBSM', 'PGSM', 'PNSM', 'POSM', 'XNSM', ...
          'DAPS', 'DOPS', 'DPPS', 'DUPS', 'OUPS', 'PAPS', 'PIPS', 'POPS', 'PQPS', 'PUPS', ...
          'PAPI', 'PIPI', 'POPI', 'PUPI', 'PAP1', 'PAP2', 'PAP3', 'POP1', 'POP2', 'POP3', ...
          'PAPA', 'PIPA', 'POPA', 'PUPA', 'PADG', 'PIDG', 'PODG', 'PUDG', 'DBCE', 'DPCE', 'DXCE', 'PNCE', 'POCE', 'XNCE', ...
           'APC',  'IPC',  'OPC',  'PPC',  'UPC',  'IPE',  'PPE', 'CHOL', ...
          'DBG1', 'DPG1', 'DXG1', 'PNG1', 'POG1', 'XNG1', 'DBG3', 'DPG3', 'DXG3', 'PNG3', 'POG3', 'XNG3', 'DBGS', 'DPGS', 'PNGS', 'POGS'};

% Lipid groups
group.name = {'PC','PE','SM','GM','GM1','GM3','GMS','CE','LPC','CHOL','DAG','PS','PI','PA','PIP','FS','MU','PU1','PU2','PU'};
% By headgroup
group.PC = {'POPX','PIPX','DPPX','DAPC', 'DOPC', 'DPPC', 'OIPC', 'OUPC', 'PAPC', 'PEPC', 'PFPC', 'PIPC', 'POPC', 'PUPC'};
group.PE = {'DAPE', 'DOPE', 'DUPE', 'OAPE', 'OIPE', 'OUPE', 'PAPE', 'PIPE', 'POPE', 'PQPE', 'PUPE'};
group.SM = {'BNSM', 'DBSM', 'DPSM', 'DXSM', 'PBSM', 'PGSM', 'PNSM', 'POSM', 'XNSM'};
group.GM = {'DBG1', 'DPG1', 'DXG1', 'PNG1', 'POG1', 'XNG1', 'DBG3', 'DPG3', 'DXG3', 'PNG3', 'POG3', 'XNG3', 'DBGS', 'DPGS', 'PNGS', 'POGS'};
group.GM1 = {'DBG1', 'DPG1', 'DXG1', 'PNG1', 'POG1', 'XNG1'};
group.GM3 = {'DBG3', 'DPG3', 'DXG3', 'PNG3', 'POG3', 'XNG3'};
group.GMS = {'DBGS', 'DPGS', 'PNGS', 'POGS'};
group.CE = {'DBCE', 'DPCE', 'DXCE', 'PNCE', 'POCE', 'XNCE'};
group.LPC = {'APC', 'IPC', 'OPC', 'PPC', 'UPC', 'IPE', 'PPE'};
group.DAG = {'PODG','PIDG', 'PADG', 'PUDG'};
group.CHOL = {'CHOL'};
group.PS = {'DAPS', 'DOPS', 'DPPS', 'DUPS', 'OUPS', 'PAPS', 'PIPS', 'POPS', 'PQPS', 'PUPS'};
group.PI = {'POPI', 'PIPI', 'PAPI', 'PUPI'};
group.PA = {'POPA', 'PIPA', 'PAPA', 'PUPA'};
group.PIP = {'PAP1', 'PAP2', 'PAP3', 'POP1', 'POP2', 'POP3'};
% By tail saturation
group.FS = {'DPPX', 'DPPC', 'DBSM', 'DPSM', 'DXSM', 'PBSM', 'DPPS', 'DBCE', 'DPCE', 'DXCE', 'PPC', 'PPE', 'DBG1', 'DPG1', 'DXG1', 'DBG3', 'DPG3', 'DXG3', 'DBGS', 'DPGS'};
group.MU = {'POPX', 'DOPC', 'POPC', 'DOPE', 'POPE', 'BNSM', 'PGSM', 'PNSM', 'POSM', 'XNSM', 'DOPS', 'POPS', 'POPI', 'POP1', 'POP2', 'POP3', 'POPA', 'PODG', 'PNCE', 'POCE', 'XNCE', 'OPC', 'PNG1', 'POG1', 'XNG1', 'PNG3', 'POG3', 'XNG3', 'PNGS', 'POGS'};
group.PU1 = {'PIPX', 'OIPC', 'OUPC', 'PAPC', 'PEPC', 'PFPC', 'PIPC', 'PUPC', 'OAPE', 'OIPE', 'OUPE', 'PAPE', 'PIPE', 'PQPE', 'PUPE', 'OUPS', 'PAPS', 'PIPS', 'PQPS', 'PUPS', 'PAPI', 'PIPI', 'PUPI', 'PAP1', 'PAP2', 'PAP3', 'PAPA', 'PIPA', 'PUPA', 'PADG', 'PIDG', 'PUDG', 'APC', 'IPC', 'UPC', 'IPE'};
group.PU2 = {'DAPC', 'DUPE', 'DAPE', 'DAPS', 'DUPS'};
group.PU = [group.PU1, group.PU2];

%% Load data
% All variables are loaded into matlab structures memsurfb (for bottom leaflet), 
% memsurft (for top leaflet), and boxsize (for the box size at a given frame)

% Loop over membrane type (APM-dep, APM-hyp, BPM-dep, BPM-hyp)
for ii = 1:length(PM_list)
    PM = PM_list{ii};
    
    % Loop over different simulation runs
    for jj = 1:length(sim_list)
        sim = sim_list{jj};
        
        % Loop over individual membranes
        for mem = 1%:4
            clear memsurfb memsurft boxsize
            
            % Loop over selected frames of the trajectory
            for i = frames+1
                if isfile([PM,'/mem',num2str(mem),'/',sim,'/','box','_',num2str(i-1),'.txt'])
                    % Load data on box size
                    boxsize{i} = textread([PM,'/mem',num2str(mem),'/',sim,'/','box','_',num2str(i-1),'.txt']);
                    
                    % Load variables extracted by MemSurfer
                    for j = 1:length(varname)
                        memsurfb.(varname{j}){i} = textread([PM,'/mem',num2str(mem),'/',sim,'/',varname{j},'b','_',num2str(i-1),'.txt']);
                        memsurft.(varname{j}){i} = textread([PM,'/mem',num2str(mem),'/',sim,'/',varname{j},'t','_',num2str(i-1),'.txt']);
                        % Correct cosine of the dipole angle for the bottom leaflet
                        if contains(varname{j},'dipl')
                            memsurfb.(varname{j}){i} = -memsurfb.(varname{j}){i};
                        end

                    end
                    % Order parameter computed for each angle, then averaged
                    memsurfb.('P2i_mean'){i} = nanmean(0.5*(3*memsurfb.('order'){i}.^2 -1),2);
                    memsurft.('P2i_mean'){i} = nanmean(0.5*(3*memsurft.('order'){i}.^2 -1),2);
                    % Order parameter computed for the average of all angles
                    memsurfb.('P2_cosTAmean'){i} = 0.5*(3*nanmean(memsurfb.('order'){i},2).^2 -1);
                    memsurft.('P2_cosTAmean'){i} = 0.5*(3*nanmean(memsurft.('order'){i},2).^2 -1);
                    % Lipid lablels  = lipid names
                    memsurfb.('label'){i} = textread([PM,'/mem',num2str(mem),'/',sim,'/','label','b','_',num2str(i-1),'.txt'],'%s');
                    memsurft.('label'){i} = textread([PM,'/mem',num2str(mem),'/',sim,'/','label','t','_',num2str(i-1),'.txt'],'%s');
                                        
                    % Load data on lipid groups
                    for j = 1:length(group.name)
                        g = group.name{j};
                        memsurfb.(g){i} = ismember(memsurfb.label{i},group.(g));
                        memsurft.(g){i} = ismember(memsurft.label{i},group.(g));
                    end
                end
            end
            % Save data into .mat file
            save([PM,'/mem',num2str(mem),'/',sim,'/','memsurf','.mat'],'memsurfb','memsurft','boxsize','group')
        end
    end
end
