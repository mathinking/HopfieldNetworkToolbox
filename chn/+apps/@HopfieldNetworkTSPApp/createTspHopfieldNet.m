function createTspHopfieldNet(app,varargin)

    app.net = [];
    C = str2double(app.parameterCEdit.String);
    
    if app.simFcnMenu.Value > 1
        simFcn = app.simFcnMenu.String{app.simFcnMenu.Value};
    else 
        simFcn = 'talavan-yanez';
    end    

    if isequal(app.tabgroup.SelectedTab,app.tabTSPLIB)
        problem = varargin{1};
        options = tsphopfieldnetOptions('Coordinates',problem.Coordinates,...
            'DistanceMatrix',problem.DistanceMatrix,'DistanceType',problem.DistanceType,'SimFcn',simFcn);
        app.net = tsphopfieldnet(problem.NumberOfCities, C, options);        
    elseif isequal(app.tabgroup.SelectedTab,app.tabPolygon) 
        N = varargin{1};
        options = tsphopfieldnetOptions('DistanceType','EUC','SimFcn',simFcn);
        app.net = tsphopfieldnet(N, C, options);
    elseif isequal(app.tabgroup.SelectedTab,app.tabOther)
        if isequal(app.tabOtherCoordsOrDist.SelectedObject,app.tabOtherCoordinates)
            coords = varargin{1};
            N = size(coords,1);
            options = tsphopfieldnetOptions('Coordinates',coords,'DistanceType','EUC','SimFcn',simFcn);
        else
            d = varargin{1};
            N = size(d,1);
            options = tsphopfieldnetOptions('DistanceMatrix',d,'DistanceType','EXPLICIT','SimFcn',simFcn);
        end
        app.net = tsphopfieldnet(N, C, options);
    end

    train(app.net);
    trainParams = getTrainParam(app.net);

    app.parameterAEdit.String = num2str(trainParams.A);
    app.parameterBEdit.String = num2str(trainParams.B);
    if strcmp(simFcn, 'talavan-yanez')
        app.parameterDEdit.String = num2str(1/trainParams.dUaux);
    else
        app.parameterDEdit.String = num2str(1/trainParams.dU);
    end
    app.parameterNEdit.String = num2str(trainParams.N);
    app.parameterNpEdit.String = num2str(trainParams.Np);
    app.cities_tauEdit.String = num2str(getCities(app.net,'Tau'));

    plot(app.net,'empty',[],app.plot);

    % Setting default for energyplot (no simulation has been carried out yet)
    if ~isempty(app.energyplot.Children)
        app.energyplot.Children.XData = [];
        app.energyplot.Children.YData = [];

        app.energyplot.XLim = [0,1];
        app.energyplot.YLim = [0,1];
        app.energyplot.XTick = 0:0.2:1;
        app.energyplot.YTick = 0:0.2:1;
        app.energyplot.YScale = 'linear';
    end

end
