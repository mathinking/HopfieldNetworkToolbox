function createFigure(app)
    
    warning('off','tsphopfieldnet:NotSimulated');
    warning('off','MATLAB:MKDIR:DirectoryExists');
    
    app.figure                  = figure('Visible','off');
    app.figure.ToolBar          = 'None';
    app.figure.MenuBar          = 'None';
    app.figure.Units            = 'normalized';
    app.figure.Position         = [0.2,0.2,0.6,0.6];
    app.figure.Interruptible    = 'off';
    app.figure.NumberTitle      = 'off';
    app.figure.Resize           = 'on';
    app.figure.HandleVisibility	= 'off';
    app.figure.Name             = 'Hopfield Net TSP solver App';
    app.figure.CloseRequestFcn  = @(src,evt)cb_closeApp(app);
    
    problems = utils.TSPLIB.problemNames(true,false);
    
    app.tabgroup = uitabgroup(app.figure);
    app.tabgroup.Position = [0.05,0.7,0.45,0.25];
    app.tabTSPLIB = uitab(app.tabgroup);
    app.tabTSPLIB.Title = 'TSPLIB';
    app.tabTSPLIB.ButtonDownFcn = @(~,~)cleanState(app);

    app.tabTSPLIBmenu = uicontrol(app.tabTSPLIB,'Style','popupmenu');
    app.tabTSPLIBmenu.Units = 'normalized';
    app.tabTSPLIBmenu.Position = [0.06,0.55,0.5,0.2];
	app.tabTSPLIBmenu.FontUnits = 'normalized';
    app.tabTSPLIBmenu.FontSize = 0.45;
    app.tabTSPLIBmenu.String = ' ';
    app.tabTSPLIBmenu.Callback = @(~,~)cb_tabTSPLIBmenu(app);
    app.tabTSPLIBmenu.String = ['Choose a TSPLIB problem...';problems];

    app.tabTSPLIBnCities = uicontrol(app.tabTSPLIB,'Style','text');
    app.tabTSPLIBnCities.Units = 'normalized';
    app.tabTSPLIBnCities.Position = [0.6,0.52,0.26,0.2];
	app.tabTSPLIBnCities.FontUnits = 'normalized';
    app.tabTSPLIBnCities.FontSize = 0.45;
    app.tabTSPLIBnCities.String = 'Number of cities:';
    
    app.tabTSPLIBnCitiesEdit = uicontrol(app.tabTSPLIB,'Style','text');
    app.tabTSPLIBnCitiesEdit.Units = 'normalized';
    app.tabTSPLIBnCitiesEdit.Position = [0.86,0.52,0.1,0.2];
	app.tabTSPLIBnCitiesEdit.FontUnits = 'normalized';
    app.tabTSPLIBnCitiesEdit.FontSize = 0.45;
    app.tabTSPLIBnCitiesEdit.String = '- - -';

    app.tabTSPLIBdistType = uicontrol(app.tabTSPLIB,'Style','text');
    app.tabTSPLIBdistType.Units = 'normalized';
    app.tabTSPLIBdistType.Position = [0.325,0.2,0.25,0.2];
	app.tabTSPLIBdistType.FontUnits = 'normalized';
    app.tabTSPLIBdistType.FontSize = 0.45;
    app.tabTSPLIBdistType.String = 'Distance Type:';
    
    app.tabTSPLIBdistTypeEdit = uicontrol(app.tabTSPLIB,'Style','text');
    app.tabTSPLIBdistTypeEdit.Units = 'normalized';
    app.tabTSPLIBdistTypeEdit.Position = [0.575,0.2,0.2,0.2];
	app.tabTSPLIBdistTypeEdit.FontUnits = 'normalized';
    app.tabTSPLIBdistTypeEdit.FontSize = 0.45;
    app.tabTSPLIBdistTypeEdit.HorizontalAlignment = 'left';  
    app.tabTSPLIBdistTypeEdit.String = '- - -';
   
    app.tabPolygon = uitab(app.tabgroup);
    app.tabPolygon.Title = 'Polygon';
    app.tabPolygon.ButtonDownFcn = @(~,~)cleanState(app);
    
    app.tabPolygonnCities = uicontrol(app.tabPolygon,'Style','text');
    app.tabPolygonnCities.Units = 'normalized';
    app.tabPolygonnCities.Position = [0.25,0.52,0.3,0.2];
    app.tabPolygonnCities.String = 'Number of cities:';
	app.tabPolygonnCities.FontUnits = 'normalized';
    app.tabPolygonnCities.FontSize = 0.45;
    
    app.tabPolygonnCitiesEdit = uicontrol(app.tabPolygon,'Style','edit');
    app.tabPolygonnCitiesEdit.Units = 'normalized';
    app.tabPolygonnCitiesEdit.Position = [0.55,0.56,0.2,0.2];
	app.tabPolygonnCitiesEdit.FontUnits = 'normalized';
    app.tabPolygonnCitiesEdit.FontSize = 0.5;
	app.tabPolygonnCitiesEdit.Callback = @(~,~)cb_tabPolygonNcities(app);
        
    app.tabPolygondistType = uicontrol(app.tabPolygon,'Style','text');
    app.tabPolygondistType.Units = 'normalized';
    app.tabPolygondistType.Position = [0.25,0.2,0.5,0.2];
	app.tabPolygondistType.FontUnits = 'normalized';
    app.tabPolygondistType.FontSize = 0.45;    
    app.tabPolygondistType.String = 'Distance Type: EUC';

    app.tabOther = uitab(app.tabgroup);
    app.tabOther.Title = 'Other TSP problem';
    app.tabOther.ButtonDownFcn = @(~,~)cleanState(app);
    
    app.tabOtherCoordsOrDist = uibuttongroup(app.tabOther);
    app.tabOtherCoordsOrDist.Units = 'normalized';   
    app.tabOtherCoordsOrDist.Position = [0.06,0.3,0.44,0.6];
    app.tabOtherCoordsOrDist.Title = 'Coordinates or Distance Matrix';
    app.tabOtherCoordsOrDist.FontSize = 9;
    app.tabOtherCoordsOrDist.FontUnits = 'normalized';
    app.tabOtherCoordsOrDist.FontAngle = 'italic';
        
    app.tabOtherCoordinates = uicontrol(app.tabOtherCoordsOrDist,'Style','radiobutton');
    app.tabOtherCoordinates.Units = 'normalized';
    app.tabOtherCoordinates.Position = [0.15,0.55,0.7,0.35];
    app.tabOtherCoordinates.String = 'Coordinates';   
    app.tabOtherCoordinates.FontUnits = 'normalized';
    app.tabOtherCoordinates.FontSize = 0.5;
    app.tabOtherCoordinates.Callback = @(~,~)cb_tabOtherCoordsDist(app);

    app.tabOtherDistance = uicontrol(app.tabOtherCoordsOrDist,'Style','radiobutton');
    app.tabOtherDistance.Units = 'normalized';
    app.tabOtherDistance.Position = [0.15,0.15,0.7,0.35];
    app.tabOtherDistance.String = 'Distance Matrix';   
    app.tabOtherDistance.FontUnits = 'normalized';
    app.tabOtherDistance.FontSize = 0.5;
	app.tabOtherDistance.Callback = @(~,~)cb_tabOtherCoordsDist(app);

    app.tabOtherdistType = uicontrol(app.tabOther,'Style','text');
    app.tabOtherdistType.Units = 'normalized';
    app.tabOtherdistType.Position = [0.06,0.15,0.44,0.1];
    app.tabOtherdistType.FontUnits = 'normalized';
    app.tabOtherdistType.FontSize = 0.85;
    app.tabOtherdistType.String = 'Distance Type: - - -';
    
    app.tabOthernCities = uicontrol(app.tabOther,'Style','text');
    app.tabOthernCities.Units = 'normalized';
    app.tabOthernCities.Position = [0.575,0.66,0.35,0.2];
	app.tabOthernCities.FontUnits = 'normalized';
    app.tabOthernCities.FontSize = 0.45;
    app.tabOthernCities.String = 'Number of cities: - - -';

    app.tabOtherChooseFromWs = uicontrol(app.tabOther,'Style','listbox');
    app.tabOtherChooseFromWs.Units = 'normalized';
    app.tabOtherChooseFromWs.Position = [0.575,0.15,0.375,0.55];
	app.tabOtherChooseFromWs.FontUnits = 'normalized';
    app.tabOtherChooseFromWs.FontSize = 0.15;    
    app.tabOtherChooseFromWs.Callback = @(~,~) cb_tabOtherChooseFromWs(app);
    app.tabOtherChooseFromWs.String = {'Choose from workspace'}; 

    % ---- Simulation ---- %
    app.simulation = uibuttongroup(app.figure);
    app.simulation.Position = [0.05,0.15,0.25,0.525];
    app.simulation.Title = 'Simulation';
    app.simulation.FontSize = 11;
    app.simulation.FontUnits = 'normalized';
    app.simulation.FontAngle = 'italic';
    app.simulation.FontWeight = 'bold';

    app.schemeMenu = uicontrol(app.simulation,'Style','popupmenu');
    app.schemeMenu.Units = 'normalized';
    app.schemeMenu.Position = [0.1,0.85,0.45,0.1];
    app.schemeMenu.FontSize = 8;
    app.schemeMenu.FontUnits = 'normalized';
    app.schemeMenu.String = ' ';
    app.schemeMenu.Callback = @(~,~)cb_scheme(app);
    app.schemeMenu.String = {'Scheme';'classic';'classic&2opt';'divide-conquer';'divide-conquer&2opt'};
    
    app.simFcnMenu = uicontrol(app.simulation,'Style','popupmenu');
    app.simFcnMenu.Units = 'normalized';
    app.simFcnMenu.Position = [0.575,0.85,0.325,0.1];
    app.simFcnMenu.FontSize = 8;
    app.simFcnMenu.FontUnits = 'normalized';
    app.simFcnMenu.String = ' ';
    app.simFcnMenu.Callback = @(~,~)cb_simFcn(app);
    app.simFcnMenu.String = {'Algorithm';'euler';'runge-kutta';'talavan-yanez'};
    
    app.ExecutionEnvironment = uibuttongroup(app.simulation);
    app.ExecutionEnvironment.Units = 'normalized';   
    app.ExecutionEnvironment.Position = [0.1,0.675,0.8,0.175];
    app.ExecutionEnvironment.Title = 'Processor';
    app.ExecutionEnvironment.FontSize = 9;
    app.ExecutionEnvironment.FontUnits = 'normalized';
    app.ExecutionEnvironment.FontAngle = 'italic';
    
    app.ExecutionEnvironmentCPU = uicontrol(app.ExecutionEnvironment,'Style','radiobutton');
	app.ExecutionEnvironmentCPU.Units = 'normalized';
    app.ExecutionEnvironmentCPU.Position = [0.15,0.35,0.3,0.35];
    app.ExecutionEnvironmentCPU.String = 'CPU';
	app.ExecutionEnvironmentCPU.FontUnits = 'normalized';
    app.ExecutionEnvironmentCPU.FontSize = 0.7;  

    app.ExecutionEnvironmentGPU = uicontrol(app.ExecutionEnvironment,'Style','radiobutton');
    app.ExecutionEnvironmentGPU.Units = 'normalized';
    app.ExecutionEnvironmentGPU.Position = [0.6,0.35,0.3,0.35];
    app.ExecutionEnvironmentGPU.String = 'GPU';   
    app.ExecutionEnvironmentGPU.FontUnits = 'normalized';
    app.ExecutionEnvironmentGPU.FontSize = 0.7;

    app.seed = uibuttongroup(app.simulation);
    app.seed.Units = 'normalized';
    app.seed.Position = [0.1,0.475,0.8,0.175];
    app.seed.Title = 'Seed';
    app.seed.FontSize = 9;
    app.seed.FontUnits = 'normalized';
    app.seed.FontAngle = 'italic';    
    
    app.seedShuffle = uicontrol(app.seed,'Style','radiobutton');
	app.seedShuffle.Units = 'normalized';
    app.seedShuffle.Position = [0.15,0.35,0.35,0.35];
    app.seedShuffle.String = 'Shuffled';
	app.seedShuffle.FontUnits = 'normalized';
    app.seedShuffle.FontSize = 0.7;  
    app.seedShuffle.Callback = @(~,~)cb_seedShuffle(app);

    app.seedFixed = uicontrol(app.seed,'Style','radiobutton');
    app.seedFixed.Units = 'normalized';
    app.seedFixed.Position = [0.6,0.35,0.3,0.35];
    app.seedFixed.String = 'Fixed';   
    app.seedFixed.FontUnits = 'normalized';
    app.seedFixed.FontSize = 0.7;
    app.seedFixed.Callback = @(~,~)cb_seedFixed(app);
    
	app.seedFixedEdit = uicontrol(app.seed,'Style','edit');
    app.seedFixedEdit.Units = 'normalized';
    app.seedFixedEdit.Position = [0.85,0.35,0.1,0.375];
    app.seedFixedEdit.FontUnits = 'normalized';
    app.seedFixedEdit.Enable = 'off';
    
    app.settings = uibuttongroup(app.simulation);
    app.settings.Units = 'normalized';
    app.settings.Position = [0.1,0.05,0.8,0.4];
    app.settings.Title = 'Settings';
    app.settings.FontSize = 9;
    app.settings.FontUnits = 'normalized';
    app.settings.FontAngle = 'italic';    
    
    app.settings_u0 = uicontrol(app.settings,'Style','text');
    app.settings_u0.Units = 'normalized';
    app.settings_u0.Position = [0.05,0.75,0.25,0.15];
    app.settings_u0.FontUnits = 'normalized';
    app.settings_u0.String = 'u0';
    app.settings_u0.HorizontalAlignment = 'right';
    
	app.settings_u0Edit = uicontrol(app.settings,'Style','edit');
    app.settings_u0Edit.Units = 'normalized';
    app.settings_u0Edit.Position = [0.3,0.76,0.225,0.15];
    app.settings_u0Edit.FontUnits = 'normalized';
    app.settings_u0Edit.Enable = 'off';
    
    app.settings_maxIter = uicontrol(app.settings,'Style','text');
    app.settings_maxIter.Units = 'normalized';
    app.settings_maxIter.Position = [0.05,0.52,0.25,0.15];
    app.settings_maxIter.FontUnits = 'normalized';
    app.settings_maxIter.String = 'MaxIter';
    app.settings_maxIter.HorizontalAlignment = 'right';
    
	app.settings_maxIterEdit = uicontrol(app.settings,'Style','edit');
    app.settings_maxIterEdit.Units = 'normalized';
    app.settings_maxIterEdit.Position = [0.3,0.53,0.225,0.15];
    app.settings_maxIterEdit.FontUnits = 'normalized';
    app.settings_maxIterEdit.Enable = 'off';
    
    app.settings_R_ITER = uicontrol(app.settings,'Style','text');
    app.settings_R_ITER.Units = 'normalized';
    app.settings_R_ITER.Position = [0.05,0.29,0.25,0.15];
    app.settings_R_ITER.FontUnits = 'normalized';
    app.settings_R_ITER.String = 'R_ITER';
    app.settings_R_ITER.HorizontalAlignment = 'right';
    
	app.settings_R_ITEREdit = uicontrol(app.settings,'Style','edit');
    app.settings_R_ITEREdit.Units = 'normalized';
    app.settings_R_ITEREdit.Position = [0.3,0.3,0.225,0.15];
    app.settings_R_ITEREdit.FontUnits = 'normalized';
    app.settings_R_ITEREdit.Enable = 'off';
    
    app.settings_dt = uicontrol(app.settings,'Style','text');
    app.settings_dt.Units = 'normalized';
    app.settings_dt.Position = [0.6,0.75,0.1,0.15];
    app.settings_dt.FontUnits = 'normalized';
    app.settings_dt.String = 'dt';
    app.settings_dt.HorizontalAlignment = 'right';
    
	app.settings_dtEdit = uicontrol(app.settings,'Style','edit');
    app.settings_dtEdit.Units = 'normalized';
    app.settings_dtEdit.Position = [0.7,0.76,0.2,0.15];
    app.settings_dtEdit.FontUnits = 'normalized';
    app.settings_dtEdit.Enable = 'off';

    app.settings_e = uicontrol(app.settings,'Style','text');
    app.settings_e.Units = 'normalized';
    app.settings_e.Position = [0.6,0.52,0.1,0.15];
    app.settings_e.FontUnits = 'normalized';
    app.settings_e.String = 'e';
    app.settings_e.HorizontalAlignment = 'right';    

	app.settings_eEdit = uicontrol(app.settings,'Style','edit');
    app.settings_eEdit.Units = 'normalized';
    app.settings_eEdit.Position = [0.7,0.53,0.2,0.15];
    app.settings_eEdit.FontUnits = 'normalized';
    app.settings_eEdit.Enable = 'off';

	app.settings_q = uicontrol(app.settings,'Style','text');
    app.settings_q.Units = 'normalized';
    app.settings_q.Position = [0.6,0.29,0.1,0.15];
    app.settings_q.FontUnits = 'normalized';
    app.settings_q.String = 'q';
    app.settings_q.HorizontalAlignment = 'right';    

	app.settings_qEdit = uicontrol(app.settings,'Style','edit');
    app.settings_qEdit.Units = 'normalized';
    app.settings_qEdit.Position = [0.7,0.3,0.2,0.15];
    app.settings_qEdit.FontUnits = 'normalized'; 
    app.settings_qEdit.Enable = 'off';
    
	app.cities_tau = uicontrol(app.settings,'Style','text');
    app.cities_tau.Units = 'normalized';
    app.cities_tau.Position = [0.05,0.06,0.25,0.15];
    app.cities_tau.FontUnits = 'normalized';
    app.cities_tau.String = 'tau';
    app.cities_tau.HorizontalAlignment = 'right';    

	app.cities_tauEdit = uicontrol(app.settings,'Style','edit');
    app.cities_tauEdit.Units = 'normalized';
    app.cities_tauEdit.Position = [0.3,0.07,0.225,0.15];
    app.cities_tauEdit.FontUnits = 'normalized';     
    app.cities_tauEdit.Enable = 'off';
    
    app.findTourPushbutton = uicontrol(app.figure,'Style','pushbutton');
    app.findTourPushbutton.Units = 'normalized';
    app.findTourPushbutton.Position = [0.125,0.075,0.1,0.05];
    app.findTourPushbutton.FontUnits = 'normalized';
    app.findTourPushbutton.String = 'Find Tour';
    app.findTourPushbutton.FontSize = 0.4;
    app.findTourPushbutton.Callback = @(~,~)cb_findTour(app);
    
    app.elapsedTime = uicontrol(app.figure,'Style','text');
    app.elapsedTime.Units = 'normalized';
    app.elapsedTime.Position = [0.05,0.025,0.25,0.0275];
    app.elapsedTime.FontUnits = 'normalized';
    app.elapsedTime.FontSize = 0.6;
    app.elapsedTime.FontAngle = 'italic';
    app.elapsedTimer = timer('Period',0.1,'TasksToExecute',inf,'ExecutionMode',...
        'fixedRate','StartFcn', @app.timeElapsed, 'TimerFcn',@app.timeElapsed);
    app.elapsedTimerCleaned = onCleanup(@()disp('delete'));% delete(app.elapsedTimer)

    % ---- Parameters ---- %
    app.parameters = uibuttongroup(app.figure);
    app.parameters.Position = [0.325,0.05,0.175,0.625];
    app.parameters.Title = 'Parameters';
    app.parameters.FontSize = 11;
	app.parameters.FontUnits = 'normalized';
    app.parameters.FontAngle = 'italic';
    app.parameters.FontWeight = 'bold';
    
    app.parameterA = uicontrol(app.parameters,'Style','text');
    app.parameterA.Units = 'normalized';   
    app.parameterA.Position = [0.15,0.825,0.2,0.075];
    app.parameterA.FontUnits = 'normalized';   
    app.parameterA.String = 'A';
    app.parameterA.FontSize = 0.5;

	app.parameterAEdit = uicontrol(app.parameters,'Style','edit');
    app.parameterAEdit.Units = 'normalized';   
    app.parameterAEdit.Position = [0.375,0.84,0.4,0.075];
    app.parameterAEdit.FontUnits = 'normalized';   
    app.parameterAEdit.FontSize = 0.5;
    app.parameterAEdit.Enable = 'off';

    app.parameterB = uicontrol(app.parameters,'Style','text');
    app.parameterB.Units = 'normalized';   
    app.parameterB.Position = [0.15,0.675,0.2,0.075];
    app.parameterB.FontUnits = 'normalized';   
    app.parameterB.String = 'B';
    app.parameterB.FontSize = 0.5;
    
    app.parameterBEdit = uicontrol(app.parameters,'Style','edit');
    app.parameterBEdit.Units = 'normalized';   
    app.parameterBEdit.Position = [0.375,0.69,0.4,0.075];
    app.parameterBEdit.FontUnits = 'normalized';   
    app.parameterBEdit.FontSize = 0.5;
    app.parameterBEdit.Enable = 'off';
    
    app.parameterC = uicontrol(app.parameters,'Style','text');
    app.parameterC.Units = 'normalized';   
    app.parameterC.Position = [0.15,0.525,0.2,0.075];
    app.parameterC.FontUnits = 'normalized';   
    app.parameterC.String = 'C';
    app.parameterC.FontSize = 0.5;

    app.parameterCEdit = uicontrol(app.parameters,'Style','edit');
    app.parameterCEdit.Units = 'normalized';   
    app.parameterCEdit.Position = [0.375,0.54,0.4,0.075];
    app.parameterCEdit.FontUnits = 'normalized';   
    app.parameterCEdit.FontSize = 0.5;
    app.parameterCEdit.Callback = @(~,~)cb_parameterCEdit(app);
    
    app.parameterD = uicontrol(app.parameters,'Style','text');
    app.parameterD.Units = 'normalized';   
    app.parameterD.Position = [0.15,0.375,0.2,0.075];
    app.parameterD.FontUnits = 'normalized';   
    app.parameterD.String = 'D';
    app.parameterD.FontSize = 0.5;

    app.parameterDEdit = uicontrol(app.parameters,'Style','edit');
    app.parameterDEdit.Units = 'normalized';   
    app.parameterDEdit.Position = [0.375,0.39,0.4,0.075];
    app.parameterDEdit.FontUnits = 'normalized';   
    app.parameterDEdit.FontSize = 0.5;
    app.parameterDEdit.Enable = 'off';
        
    app.parameterN = uicontrol(app.parameters,'Style','text');
    app.parameterN.Units = 'normalized';   
    app.parameterN.Position = [0.15,0.225,0.2,0.075];
    app.parameterN.FontUnits = 'normalized';   
    app.parameterN.String = 'N';
    app.parameterN.FontSize = 0.5;

    app.parameterNEdit = uicontrol(app.parameters,'Style','edit');
    app.parameterNEdit.Units = 'normalized';   
    app.parameterNEdit.Position = [0.375,0.24,0.4,0.075];
    app.parameterNEdit.FontUnits = 'normalized';   
    app.parameterNEdit.FontSize = 0.5;
    app.parameterNEdit.Enable = 'off';
    
    app.parameterNp = uicontrol(app.parameters,'Style','text');
    app.parameterNp.Units = 'normalized';   
    app.parameterNp.Position = [0.15,0.075,0.2,0.075];
    app.parameterNp.FontUnits = 'normalized';   
    app.parameterNp.String = 'N''';
    app.parameterNp.FontSize = 0.5;

    app.parameterNpEdit = uicontrol(app.parameters,'Style','edit');
    app.parameterNpEdit.Units = 'normalized';   
    app.parameterNpEdit.Position = [0.375,0.09,0.4,0.075];
    app.parameterNpEdit.FontUnits = 'normalized';   
    app.parameterNpEdit.FontSize = 0.5;  
    app.parameterNpEdit.Enable = 'off';
    
    app.plot = axes('Parent',app.figure);
    app.plot.Position = [0.55,0.4,0.4,0.5];
    app.plot.Box = 'on';
    app.plot.Title.String = 'Problem Coordinates';

    app.energyplot = axes('Parent',app.figure);
    app.energyplot.Position = [0.55,0.05,0.4,0.25];    
    app.energyplot.Box = 'on';
    app.energyplot.Title.String = 'Energy Function';
    
end
