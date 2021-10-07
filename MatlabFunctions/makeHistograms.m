function [H1,H2,KLD] = makeHistograms(porated,nonporated,memids,varname,tpbt,EDGES,nbins,doPlot)
    
for mem = memids
    if doPlot
        figure;
        nrows = ceil(length(varname)/3);
        set(gcf,'position',[69 24 1375 250*nrows+150])
    end
    
    for j = 1:length(varname)

        %------------------------ BINNING ---------------------------------
        % Predefined binning 
        if numel(EDGES) > 1
            edges = EDGES{j};
        elseif numel(EDGES) == 1
            edges = EDGES;
        
        % Automatic definition of the binning
        else
            % Compute probability density distribution using kernel smoothing
            [pd1,x1] = ksdensity(porated{mem}.(tpbt).(varname{j}),'NumPoints',200);
            [pd2,x2] = ksdensity(nonporated{mem}.(tpbt).(varname{j}),'NumPoints',200);
            % Determine the minimum and maximum value of the data
            minmin = nanmin(nanmin(porated{mem}.(tpbt).(varname{j})),nanmin(nonporated{mem}.(tpbt).(varname{j})));
            maxmax = nanmax(nanmax(porated{mem}.(tpbt).(varname{j})),nanmax(nonporated{mem}.(tpbt).(varname{j})));
            % Determine the interval that covers both of the distributions; we
            % look at when the distribution falls below 1% of its peak value
            idx1a = find(pd1>=0.01*max(pd1),1,'First');
            idx1b = find(pd1>=0.01*max(pd1),1,'Last');
            idx2a = find(pd2>=0.01*max(pd2),1,'First');
            idx2b = find(pd2>=0.01*max(pd2),1,'Last');
            xmin10 = min(x1(idx1a),x2(idx2a));
            xmax10 = max(x1(idx1b),x2(idx2b));
            dist = xmax10 - xmin10;
            edges = linspace(max(minmin,xmin10-0.05*dist),min(maxmax,xmax10+0.05*dist),nbins);
            % For lipid densities for lipids that don't exist in a given leaflet
            if sum(porated{mem}.(tpbt).(varname{j})) == 0
                edges = linspace(0,1,nbins);
                doSubplot = 0;
            else
                doSubplot = 1;
            end
        end
        
        %--------------------- HISTOGRAMS ---------------------------------
        % Compute histograms on the interval "edges"
        h10 = hist(porated{mem}.(tpbt).(varname{j}),edges);
        h20 = hist(nonporated{mem}.(tpbt).(varname{j}),edges);
        % Normalized
        h1 = h10./sum(h10)./diff(edges(1:2));
        h2 = h20./sum(h20)./diff(edges(1:2));
        
        % Compute probability density distribution using kernel smoothing
        edges_pdf = linspace(edges(1),edges(end),100);
        if strcmp(varname{j},'GM') || strcmp(varname{j},'CE') || strcmp(varname{j},'LPC') || strcmp(varname{j},'DAG') || ...
                strcmp(varname{j},'PA') || strcmp(varname{j},'PIP') || strcmp(varname{j},'charg')
            [pd2,x2,U2] = ksdensity(nonporated{mem}.(tpbt).(varname{j}),edges_pdf,'BandWidth',diff(edges(1:2))/3);
            [pd1,x1,U1] = ksdensity(porated{mem}.(tpbt).(varname{j}),edges_pdf,'BandWidth',U2);
        else
            [pd2,x2,U2] = ksdensity(nonporated{mem}.(tpbt).(varname{j}),edges_pdf,'BandWidth',diff(edges(1:2)));
            [pd1,x1,U1] = ksdensity(porated{mem}.(tpbt).(varname{j}),edges_pdf,'BandWidth',U2);            
        end
        
        % -------------- Distances between histograms ---------------------
        % Take into acount only parts of the pdf which are above 1e-5 to avoid problems when computing logarithms
        idx1 = find(pd1>1e-5);
        idx2 = find(pd2>1e-5);
        pd1bin = pd1(intersect(idx1,idx2))*diff(edges_pdf(1:2));
        pd2bin = pd2(intersect(idx1,idx2))*diff(edges_pdf(1:2));
        % Symmetric kld
        kld = sum(0.5*pd1bin.*log(pd1bin./pd2bin) + 0.5*pd2bin.*log(pd2bin./pd1bin)).*sign(median(porated{mem}.(tpbt).(varname{j}))-median(nonporated{mem}.(tpbt).(varname{j})));
       
        if doPlot
            subplot(nrows,3,j); hold on; box on
            if doSubplot
                bar(edges,h1,0.95,'LineStyle','none')
                bar(edges,h2,0.95,'FaceColor','none','LineWidth',2)
                plot(x1,pd1,'LineWidth',2)
                plot(x2,pd2,'LineWidth',2)
            end
            title(strrep(varname{j},'_','\_'))
            xlim([edges(1),edges(end)])
        end
        
        % ----------------------- OUTPUTs ---------------------------------
        H1{mem}.(varname{j})  = h1;
        H2{mem}.(varname{j})  = h2;
        KLD{mem}.(varname{j}) = kld;
    end
end

