classdef tsphopfieldnetApp < handle
    
    properties (Access = public)
        figure;
        net;
        
        tabgroup;
        
        tabTSPLIB;
        tabTSPLIBmenu;
        tabTSPLIBnCities;
        tabTSPLIBnCitiesEdit;
        tabTSPLIBdistType;
        tabTSPLIBdistTypeEdit;
        
        tabPolygon;
        tabPolygonnCities;
        tabPolygonnCitiesEdit;
        tabPolygondistType;
        
        tabOther;
        tabOtherCoordsOrDist;
        tabOtherdistType;
        tabOthernCities;
        tabOtherChooseFromWs;
        tabOtherCoordinates;
        tabOtherDistance;
        
        simulation;      
        simFcnMenu;
        
        hwResources;
        hwResourcesCPU;
        hwResourcesGPU;
        
        seed;
        seedShuffle;
        seedFixed;
        seedFixedEdit;
        
        settings;
        settings_u0;
        settings_u0Edit;
        settings_maxIter;
        settings_maxIterEdit;
        settings_R_ITER;
        settings_R_ITEREdit;
        settings_dt;
        settings_dtEdit;
        settings_e;
        settings_eEdit;
        settings_q;
        settings_qEdit;
        cities_tau;
        cities_tauEdit;
        
        findTourPushbutton;
        
        elapsedTime;
        elapsedTimeSecs;
        elapsedTimer;
        elapsedTimerCleaned;
        
        parameters;
        parameterA;
        parameterB;
        parameterC;
        parameterD;
        parameterN;
        parameterNp;
        parameterAEdit;
        parameterBEdit;
        parameterCEdit;
        parameterDEdit;
        parameterNEdit;
        parameterNpEdit;
        
        plot;
        energyplot;
    end
    
    methods 
        function app = tsphopfieldnetApp()
            createFigure(app);
            defaultSettings(app);
            drawnow nocallbacks;
            app.figure.Visible = 'on';
        end
    end
       
    methods (Access = private)
        createFigure(app);
    	defaultSettings(app);
        cb_tabTSPLIBmenu(app);
        cb_tabPolygonNcities(app);
        cb_tabOtherCoordsDist(app);
        cb_tabOtherChooseFromWs(app);
        cb_seedShuffle(app);
        cb_seedFixed(app);
        cb_parameterCEdit(app);
        cb_findTour(app);
        cb_simFcn(app);
        timeElapsed(app,obj,event);
        cleanState(app);
        createTspHopfieldNet(app,varargin);
        updateTspHopfieldNet(app,C);
    end
end
