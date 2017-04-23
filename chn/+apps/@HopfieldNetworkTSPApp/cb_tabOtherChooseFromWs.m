function cb_tabOtherChooseFromWs(app)
    
    workspaceNames = app.tabOtherChooseFromWs.String;
    if app.tabOtherChooseFromWs.Value == 1
        app.tabOthernCities.String = 'Number of cities: - - -';
    else
        coordsOrDistance = evalin('base',workspaceNames{app.tabOtherChooseFromWs.Value});
        app.tabOthernCities.String = ['Number of cities: ', num2str(size(coordsOrDistance,1))];
    end
    
    % Obtain information and update app values for training and simulation
    if app.tabOtherChooseFromWs.Value > 1
        createTspHopfieldNet(app,coordsOrDistance);     
    else
        state = app.tabOtherCoordsOrDist.SelectedObject;
        defaultSettings(app);
        cleanState(app);
        app.tabOtherCoordsOrDist.SelectedObject = state;
        cb_tabOtherCoordsDist(app)
    end    
end
