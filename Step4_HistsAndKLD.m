%% Make histograms and determine the distance between histograms by KL divergence
addpath('../MatlabFunctions')
set(0,'defaultAxesFontSize',14)

%%
PM = 'APM-dep';
rr = 15; 

%% variables A
clear dist
load([PM,'/forML_equil10nsBeforeEField','_rr',num2str(rr),'.mat'])
varnameA = {'parea','thick','mcurv','dipln','charg','P2i_mean','P2_cosTAmean'};

% Plot histograms for the top leaflet
[~,~,dist.KLD1] = makeHistograms(poratedA,nonporatedA,5,varnameA,'tp',[],31,1);
% Plot histograms for the bottom leaflet
[~,~,dist.KLD2] = makeHistograms(poratedA,nonporatedA,5,varnameA,'bt',[],31,1);
% Plot histograms for the mean of both leaflets 
[~,~,dist.KLD3] = makeHistograms(poratedA,nonporatedA,5,varnameA,'mean',[],31,1);
% Plot histograms for the difference between the top and bottom leaflet
[~,~,dist.KLD4] = makeHistograms(poratedA,nonporatedA,5,varnameA,'diff',[],31,1);

% Plot distances between histograms
figure;
for ii = 1:4    
    subplot(1,4,ii); hold on; box on
    set(gca,'XTick',1:length(varnameA),'XTickLabel',strrep(varnameA,'_',' '))
    ylim([-1 1]); ylabel('KL divergence')
    for j = 1:length(varnameA)        
        bar(j,dist.(['KLD',num2str(ii)]){5}.(varnameA{j}))
    end
    switch ii
        case 1
            title('top leaflet')
        case 2
            title('bottom leaflet')
        case 3
            title('mean of leaflets')
        case 4
            title('diff of leaflets')
    end
end

%% variables B
clear dist
load([PM,'/forML_equil10nsBeforeEField','_rr',num2str(rr),'.mat'])
varnameB = {'PC','PE','SM','GM','CE','LPC','CHOL','DAG','PS','PI','PA','PIP','FS','MU','PU'};

% Plot histograms for the top leaflet
[~,~,dist.KLD1] = makeHistograms(poratedB,nonporatedB,5,varnameB,'tp',[],31,1);
% Plot histograms for the bottom leaflet
[~,~,dist.KLD2] = makeHistograms(poratedB,nonporatedB,5,varnameB,'bt',[],31,1);
% Plot histograms for the mean of both leaflets 
[~,~,dist.KLD3] = makeHistograms(poratedB,nonporatedB,5,varnameB,'mean',[],31,1);
% Plot histograms for the difference between the top and bottom leaflet
[~,~,dist.KLD4] = makeHistograms(poratedB,nonporatedB,5,varnameB,'diff',[],31,1);

% Plot distances between histograms
figure;
for ii = 1:4    
    subplot(1,4,ii); hold on; box on
    set(gca,'XTick',1:length(varnameB),'XTickLabel',strrep(varnameB,'_',' '))
    ylim([-1 1]); ylabel('KL divergence')
    for j = 1:length(varnameB)        
        bar(j,dist.(['KLD',num2str(ii)]){5}.(varnameB{j}))
    end
    switch ii
        case 1
            title('top leaflet')
        case 2
            title('bottom leaflet')
        case 3
            title('mean of leaflets')
        case 4
            title('diff of leaflets')
    end
end