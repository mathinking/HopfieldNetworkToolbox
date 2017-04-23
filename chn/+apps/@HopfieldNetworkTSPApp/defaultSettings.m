function defaultSettings(app)
% Set Default Settings for App

    app.tabTSPLIBmenu.Value = 1;
    app.tabTSPLIBnCitiesEdit.String  = '- - -';
    app.tabTSPLIBdistTypeEdit.String = '- - -';

    app.tabPolygonnCitiesEdit.String = '';

    app.tabOtherCoordsOrDist.SelectedObject = app.tabOtherCoordinates;
    app.tabOthernCities.String = 'Number of cities: - - -';
    app.tabOtherChooseFromWs.Value = 1;
    
    % Simulation
    app.simFcnMenu.Value = 1;
    app.hwResources.SelectedObject = app.hwResourcesCPU;

    app.seed.SelectedObject = app.seedShuffle;
    app.seedFixedEdit.String = '3';

    app.settings_u0Edit.String = '0.3';
    app.settings_maxIterEdit.String = '2000';
    app.settings_R_ITEREdit.String = '20';
    app.settings_dtEdit.String = '0.01';
    app.settings_eEdit.String = '13';
    app.settings_qEdit.String = '0.8';

    app.cities_tauEdit.String = '';

    % Parameters
    app.parameterAEdit.String = '';
    app.parameterBEdit.String = '';
    app.parameterCEdit.String = '1e-5';
    app.parameterDEdit.String = '';
    app.parameterNEdit.String = '';
    app.parameterNpEdit.String = '';

end
