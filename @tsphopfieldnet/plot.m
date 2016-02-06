function h = plot(net, type, chains, ax, varargin)

	mytext = net.cities.names;
        
    if nargin == 1;
        if isempty(net.results.tourLength)
            type = 'empty';
        else
            type = 'total';
        end
        chains = [];
    end
    
    if isempty(varargin)
        myCitiesColor = [1,0,0];
        myCitiesTextColor = [0,0,0.5];
        myInsideColor = [1,0,0];
        myTitle = [];
    else
        myCitiesColor = varargin{1};
        myCitiesTextColor = varargin{2};        
        myInsideColor = varargin{3}; %[0.8,0.8,0.8];
        myTitle = varargin{4};
    end

    if nargin < 3 || isempty(ax)
        figure;
        ax = axes;
    end
        
    coords = net.cities.coords;
    if isempty(coords) % No data coords available. Just visualizing order
        warning('tsphopfieldnet:NoDataCoords','No data coords available. Visualizing order instead.');
        coords = fakeCoords(net);          
    end

    % 3 Possible plots: Empty, Total, Partial

    if strcmp(type, 'empty') % Simulation has not taken place yet. Ploting data points.
%         hold(ax,'on');
        h = plot(ax,coords(:,1),coords(:,2),...
            'color','r','marker','o','markersize',4,...
            'markerfacecolor','r',...
            'markeredgecolor',myCitiesColor,'linestyle','none');
        text(coords(:,1),coords(:,2),mytext,'fontsize',10,'color',myCitiesTextColor,'margin',20,'Clipping','on','Parent',ax);
        title(ax,'\bf Problem coordinates');
        warning('tsphopfieldnet:NotSimulated','Simulation has not taken place yet. Use ''sim(net)'' to simulate your network.');
        
    else
        if net.results.validPath && (strcmp(type,'total') || strcmp(type,'phase2')) % Total plot or phase2
            dispText = strrep(cellstr([char(mytext(net.results.visitOrder)'),...
                repmat('_{(',length(mytext),1),num2str((1:net.trainParam.N)'),repmat(')}',...
                length(mytext),1)]),' ','');

            h = plot(ax,[coords(net.results.visitOrder,1);coords(net.results.visitOrder(1),1)],...
                [coords(net.results.visitOrder,2);coords(net.results.visitOrder(1),2)],...
                'color',myCitiesColor,'marker','o','markersize',4,'markerfacecolor',myCitiesColor,...
                'markeredgecolor',myCitiesColor,'linestyle','-');

            hold(ax,'on');
            if strcmp(type,'phase2')
                for c = 1:length(chains)
                    h = plot(coords(chains{c} ,1), coords(chains{c} ,2),...
                        'color',myInsideColor,'marker','o','markersize',4,'markerfacecolor',myInsideColor,...
                        'markeredgecolor',myInsideColor,'linestyle','-');                    
                end               
            end
            text(coords(net.results.visitOrder,1),coords(net.results.visitOrder,2),...
                dispText,'fontsize',10,'color',myCitiesTextColor,'margin',20,'Clipping','on','Parent',ax);
            
            hold(ax,'off');
            
        elseif net.results.validPath && strcmp(type,'phase1') % Partial plot
            
            % Plot single points
            singlePoints = setxor(1:net.trainParam.N, [chains{:}]);
            h = plot(coords(singlePoints,1), coords(singlePoints,2),...
                'color',myCitiesColor,'marker','o','markersize',4,'markerfacecolor',myCitiesColor,...
                'markeredgecolor',myCitiesColor,'linestyle','none');
            hold on;
            for c = 1:length(chains)
                % Plot Extreme Chain Points
                plot(coords(chains{c}([1,end]) ,1), coords(chains{c}([1,end]) ,2),...
                    'color',myCitiesColor,'marker','o','markersize',4,'markerfacecolor',myCitiesColor,...
                    'markeredgecolor',myCitiesColor,'linestyle','-');
                % Plot Chain
                h = plot(coords(chains{c} ,1), coords(chains{c} ,2),...
                    'color',myInsideColor,'marker','o','markersize',4,'markerfacecolor',myInsideColor,...
                    'markeredgecolor',myInsideColor,'linestyle','-');
            end

            text(coords(:,1),coords(:,2),mytext,...
                'fontsize',10,'color',myCitiesTextColor,'margin',20,'Clipping','on');
            hold off;
            
        else
            error('tsphopfieldnet:UnvalidPath','Unvalid TSP Path');
        end
        if ~isempty(myTitle)
            title(ax,['\bf', myTitle]);
        else
            title(ax,['\bf Tour length: ', num2str(net.results.tourLength)]);
        end
        
    end
    adjustPlot(h,coords);
end

function adjustPlot(h,coords)
    xMin = min(coords(:,1)); yMin = min(coords(:,2)); 
    xMax = max(coords(:,1)); yMax = max(coords(:,2));
    xRange = max(0.5,xMax - xMin); yRange = max(0.5,yMax - yMin);
    set(get(h,'Parent'),'XLim',[xMin-xRange/10,xMax+xRange/10],'YLim', [yMin-yRange/10,yMax+yRange/10])
end

function coords = fakeCoords(net)
    coords = zeros(net.trainParam.N,2);
    coords(:,1) = 1:net.trainParam.N;
    coords(1:2:end,2) =  0.005;
    coords(2:2:end,2) = -0.005;
end
