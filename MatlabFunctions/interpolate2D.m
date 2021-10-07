function [outb,outt] = interpolate2D(boxsize,datab,datat,varname,frames,X,Y,interpolation,rr,doPlot)
% Determine the value at given (x,y) location by taking account the values
% in the surrounding vertices, either by Matlab's scatteredInterpolant or
% by using a Gaussian smoothing kernel with size rr

Xp = X(:);
Yp = Y(:);

for ii = 1:length(frames)
    i = frames(ii)+1;
    
    % Box size
    boxx = boxsize{i}(2,1);
    boxy = boxsize{i}(2,2);
    % (x,y) coordinates of vertices in the top and bottom leaflet
    xt = datat.verti{i}(:,1);
    yt = datat.verti{i}(:,2);
    xb = datab.verti{i}(:,1);
    yb = datab.verti{i}(:,2);
    % (x,y) coordinates of vertices and their pbc images
    xt_pbc = [xt;       xt-boxx; xt;      xt+boxx; xt-boxx; xt+boxx; xt-boxx; xt;      xt+boxx];
    yt_pbc = [yt;       yt-boxy; yt-boxy; yt-boxy; yt;      yt;      yt+boxy; yt+boxy; yt+boxy];
    xb_pbc = [xb;       xb-boxx; xb;      xb+boxx; xb-boxx; xb+boxx; xb-boxx; xb;      xb+boxx];
    yb_pbc = [yb;       yb-boxy; yb-boxy; yb-boxy; yb;      yb;      yb+boxy; yb+boxy; yb+boxy];
    
    % Data in vertices
    tmpb = datab.(varname){i};
    tmpb_pbc= repmat(tmpb,9,1);
    tmpt = datat.(varname){i};
    tmpt_pbc= repmat(tmpt,9,1);
    
    % If coordinates are in relative units, then multiply with box size
    if max(Xp) < 2
        xp = Xp*boxx;
        yp = Yp*boxy;
    else
        xp = Xp;
        yp = Yp;
    end
    
    if strcmp(interpolation,'scatteredInterp')
        F = scatteredInterpolant(xb_pbc,yb_pbc,tmpb_pbc,'natural','none');
        outb(:,ii) = F(xp,yp);
        F = scatteredInterpolant(xt_pbc,yt_pbc,tmpt_pbc,'natural','none');
        outt(:,ii) = F(xp,yp);
        
    elseif strcmp(interpolation,'gaussKernel')
        for j = 1:length(xp)
            idx = find(((xb_pbc-xp(j)).^2 + (yb_pbc-yp(j)).^2) < (3*rr)^2);
            weights = exp(-(0.5*(xb_pbc(idx)-xp(j)).^2/rr.^2 + 0.5*(yb_pbc(idx)-yp(j)).^2/rr.^2));
            outb(j,ii) = nansum(tmpb_pbc(idx).*weights)./nansum(weights);
            idx = find(((xt_pbc-xp(j)).^2 + (yt_pbc-yp(j)).^2) < (3*rr)^2);
            weights = exp(-(0.5*(xt_pbc(idx)-xp(j)).^2/rr.^2 + 0.5*(yt_pbc(idx)-yp(j)).^2/rr.^2));
            outt(j,ii) = nansum(tmpt_pbc(idx).*weights)./nansum(weights);
        end        
    end
    
    OUTb(:,:,ii) = reshape(outb(:,ii),size(X,1),size(X,2));
    OUTt(:,:,ii) = reshape(outt(:,ii),size(X,1),size(X,2));
    
end

if doPlot
    figure; hold on; box on
    subplot(1,3,1); set(gca,'position',[0.05,0.05,0.27,0.9])
    surfc(X,Y,squeeze(nanmean(OUTb,3)),'EdgeColor','none')
    view(0,90); axis equal; colorbar; %caxis(clim);
    title([strrep(varname,'_','\_'), ' bottom'])
    subplot(1,3,2); set(gca,'position',[0.38,0.05,0.27,0.9])
    surfc(X,Y,squeeze(nanmean(OUTt,3)),'EdgeColor','none')
    view(0,90); axis equal; colorbar; %caxis(clim);
    title([strrep(varname,'_','\_'), ' top'])
    subplot(1,3,3); set(gca,'position',[0.70,0.05,0.27,0.9])
    surfc(X,Y,0.5*(squeeze(nanmean(OUTb,3)) + squeeze(nanmean(OUTt,3))),'EdgeColor','none')
    view(0,90); axis equal; colorbar; %caxis(clim);
    title('average top, bottom')
end