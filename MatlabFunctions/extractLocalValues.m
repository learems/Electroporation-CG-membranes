function [porated, nonporated] = extractLocalValues(boxsize,memsurfb,memsurft,varname1,frames,poreloci,noporeloci,interpolation,rr)
% Extract the values of selected variables from porated and nonporated locations

% Porated and nonporated locations
if numel(poreloci)
    xp = poreloci(:,1);
    yp = poreloci(:,2);
end
if numel(noporeloci)
    xnp = noporeloci(:,1);
    ynp = noporeloci(:,2);
end

for j = 1:length(varname1)
    % Prepare struct; we will output the values from the top leaflet (tp),
    % values from the bottom leaflet (bt), sum of the values from both
    % leaflets (sum) and difference between the values in the top and
    % bottom leaflet. 
    porated.tp.(varname1{j}) = [];
    porated.bt.(varname1{j}) = [];
    porated.sum.(varname1{j}) = [];
    porated.diff.(varname1{j}) = [];
    nonporated.tp.(varname1{j}) = [];
    nonporated.bt.(varname1{j}) = [];
    nonporated.sum.(varname1{j}) = [];
    nonporated.diff.(varname1{j}) = [];
    
    if numel(poreloci)
        [outb,outt] = interpolate2D(boxsize,memsurfb,memsurft,varname1{j},frames,xp,yp,interpolation,rr,0);
        porated.bt.(varname1{j}) = [porated.bt.(varname1{j}); outb(:)];
        porated.tp.(varname1{j}) = [porated.tp.(varname1{j}); outt(:)];
        porated.sum.(varname1{j}) = [porated.sum.(varname1{j}); nansum([outb(:), outt(:)],2)];
        porated.diff.(varname1{j}) = [porated.diff.(varname1{j}); outt(:)-outb(:)];
        porated.mean.(varname1{j}) = 0.5*porated.sum.(varname1{j});
    end
    if numel(noporeloci)
        [outb,outt] = interpolate2D(boxsize,memsurfb,memsurft,varname1{j},frames,xnp,ynp,interpolation,rr,0);
        nonporated.bt.(varname1{j}) = [nonporated.bt.(varname1{j}); outb(:)];
        nonporated.tp.(varname1{j}) = [nonporated.tp.(varname1{j}); outt(:)];
        nonporated.sum.(varname1{j}) = [nonporated.sum.(varname1{j}); nansum([outb(:), outt(:)],2)];
        nonporated.diff.(varname1{j}) = [nonporated.diff.(varname1{j}); outt(:)-outb(:)];
        nonporated.mean.(varname1{j}) = 0.5*nonporated.sum.(varname1{j});
    end
end
    

