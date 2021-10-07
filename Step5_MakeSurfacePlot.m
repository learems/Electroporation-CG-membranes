%% Make a surface plot
addpath('../MatlabFunctions')
set(0,'defaultAxesFontSize',14)

%% Make a surface plot for lipid densities
xgv = 0:0.02:1;
ygv = 0:0.02:1;
[X,Y] = meshgrid(xgv,ygv);
rr = 15;

PM = 'APM-dep';

for mem = 1
    load([PM,'/mem',num2str(mem),'/equil10nsBeforeEField/','memsurf.mat'])
    [outb,outt] = interpolate2D(boxsize,memsurfb,memsurft,'PU',0,X,Y,'gaussKernel',rr,1);
    % Overlay pore locations
%     load(['PoreLocations/',PM,'.mat'])
%     hold on;
%     plot3(poreloci_rel{mem}(:,1),poreloci_rel{mem}(:,2),100*ones(size(poreloci_rel{mem}(:,1))),'wo','MarkerSize',10)
end
