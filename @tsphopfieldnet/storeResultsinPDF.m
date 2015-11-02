function storeResultsinPDF(results,simulations,N_or_Problem,CF,resultsPath)
    figure; ax = gca;
    set(gcf, 'Visible', 'Off');
    plot(ax, simulations, results.perfRatio, 'Marker', 'x', ...
        'Color', [0,0.75,0], 'MarkerSize', 10, 'LineStyle', 'none', ...
        'LineWidth', 2);
    hold('on');
    set(ax,'XLim',[simulations(1)-1,simulations(end)+1]);

    if any(~isnan(results.perfRatio))
        set(ax,'YLim',[min(results.perfRatio)-0.1,max(results.perfRatio)+0.1]);
    end

    % Plot non-feasible solutions
    nonFeasible = results.exitFlag == -1;
    if any(nonFeasible)
        nonSol = find(nonFeasible);
        plot(nonSol, repmat(min(results.perfRatio)-0.1, 1, ...
            nnz(nonFeasible)), 'Marker', 'x', 'Color', [0.75,0,0], ...
            'LineStyle', 'none', 'MarkerSize', 10, 'LineWidth', 2)
    end

    % Plot non-feasible solutions (that have reached max iterations)
    nonFeasibleMaxIterations = results.exitFlag == 0;
    if any(nonFeasibleMaxIterations)
        nonSol = find(nonFeasibleMaxIterations);
        plot(nonSol, repmat(min(results.perfRatio)-0.1, 1, ...
            nnz(nonFeasibleMaxIterations)), 'Marker', 'x', ...
            'Color', [1,0.4,0], 'LineStyle', 'none', 'MarkerSize', 10, ...
            'LineWidth', 2)
    end

    % Labeling graph
    titleStr = 'Hopfield Network Results: ';
    if isa(N_or_Problem,'tsplib')
        problem = N_or_Problem;
        titleStr = [titleStr, 'TSPLIB - ', problem.name];
        N = problem.nCities;
    else
        N = N_or_Problem;
        titleStr = [titleStr, 'Cities are vertices in Polygon'];
    end
    
    C = CF(1);
    title([titleStr, ' N = ', num2str(N), '. C = ',num2str(C),'.']);
    if length(CF) == 2
        F = CF(2);
        title([titleStr, ', N = ', num2str(N), '. C = ',num2str(C), '. F = ',num2str(F),'.']);  
    end
    
    xlabel('Simulation No.');
    ylabel('Hopfield Network Performance Ratio');

    % Plotting Simulation Statistics
    plot(get(ax,'XLim'), [nanmax(results.perfRatio), ...
                          nanmax(results.perfRatio)], 'r:')
    plot(get(ax,'XLim'), [nanmin(results.perfRatio), ...
                          nanmin(results.perfRatio)], 'b:')
    plot(get(ax,'XLim'), [nanmean(results.perfRatio), ...
                          nanmean(results.perfRatio)], 'k-.')
    plot(get(ax,'XLim'), [mode(round(results.perfRatio,2)), ...
                          mode(round(results.perfRatio,2))], 'k:')

    feasibility = {'Feasible Solution'};

    if ~all(results.exitFlag == 1)
        if any(nonFeasible)
            feasibility = {feasibility{:},'Non-Feasible Solution'};
            if any(nonFeasibleMaxIterations)
                feasibility = {feasibility{:},...
                    'Non-Feasible Solution - Max. Iterations Reached'};
            end
        elseif any(nonFeasibleMaxIterations)
            feasibility = {feasibility{:},...
                'Non-Feasible Solution - Max. Iterations Reached'};
        end
    end   

    figurePos = get(gcf,'Position');
    set(gcf,'Position',...
        [figurePos(1),figurePos(2)*0.6,figurePos(3)*1.40,figurePos(4)*1.40]);
    axesPos = get(get(gcf,'Children'),'Position');
    set(get(gcf,'Children'),'Position',...
        [axesPos(1),axesPos(2)+(axesPos(4)*0.3),axesPos(3),axesPos(4)*0.7])

    myLegend = legend(feasibility{:},...
        ['max: ',num2str(nanmax(results.perfRatio))], ...
        ['min: ',num2str(nanmin(results.perfRatio))], ...
        ['mean: ',num2str(nanmean(results.perfRatio))], ...
        ['mode: ',num2str(mode(round(results.perfRatio,2)))]);%, ...
    %     'Location', 'SouthOutside');
    set(myLegend, 'Position', ...
        [(1-axesPos(3)*0.5)/2,axesPos(2)*0.50,axesPos(3)*0.5,axesPos(4)*0.2])

    % Saving visualization
    if isa(N_or_Problem,'tsplib')
        if length(CF) == 2
            visualizationFile = ['TSPLIB_',problem.name,'_C=',num2str(C),'_F=',num2str(F),'.pdf'];
        else
            visualizationFile = ['TSPLIB_',problem.name,'_C=',num2str(C),'.pdf'];
        end
    else
        if length(CF) == 2
            visualizationFile = ['Polygon_N=',problem.name,'_C=',num2str(C),'_F=',num2str(F),'.pdf'];
        else
            visualizationFile = ['Polygon_N=',problem.name,'_C=',num2str(C),'.pdf'];
        end
    end

    % set(gcf, 'Position', get(0,'Screensize'));
%     set(gcf,'PaperPositionMode','auto')
    print(gcf,'-dpdf','-r500',fullfile(resultsPath,visualizationFile))
end
