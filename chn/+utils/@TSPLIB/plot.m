function plot(problem)

    [optimTourLength, optimTour] = findOptimumTour(problem);
    optimTour = optimTour(1:end-1);
    if ~isnan(optimTourLength)
        mytext = tsphopfieldnet.cityTextGeneration(problem.nCities);
        dispText = strrep(cellstr(...
            [char(mytext(optimTour',:)),...
            repmat('_{(',length(mytext),1),...
            num2str((1:problem.nCities)'),repmat(')}',...
            length(mytext),1)]),' ','');

            figure;
            h = plot([problem.coords(optimTour,1);...
                problem.coords(optimTour(1),1)],...
                [problem.coords(optimTour,2);...
                problem.coords(optimTour(1),2)],...
                'color',[0,0.5,0],'marker','o','markersize',4,...
                'markerfacecolor',[0,1,0],...
                'markeredgecolor',[0,0.5,0],...
                'linestyle','-');

            xMin = min(problem.coords(:,1)); yMin = min(problem.coords(:,2)); 
            xMax = max(problem.coords(:,1)); yMax = max(problem.coords(:,2));
            xRange = xMax - xMin; yRange = yMax - yMin;
            set(get(h,'Parent'),...
                'XLim',[xMin-xRange/10,xMax+xRange/10],...
                'YLim', [yMin-yRange/10,yMax+yRange/10])
            hold on;
            text(problem.coords(optimTour,1),...
                problem.coords(optimTour,2),...
                dispText,'fontsize',10,...
                'color',[0.5,0.5,0.5],'margin',20,'Clipping','on');
            title(['\bf Optimum our length: ',...
                num2str(optimTourLength)]);
            hold off;
    else
        fprintf(['Optimum tour file not available for ',...
            'problem %s.\n'],problem.name);
    end
    
end