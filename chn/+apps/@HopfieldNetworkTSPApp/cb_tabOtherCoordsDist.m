function cb_tabOtherCoordsDist(app)
    w = evalin('base','whos');

    % Look for coordinates or distance matrices
    coordinates = cell(size(w));
    distances = cell(size(w));
    for i = 1:size(w,1)
        if strcmp((w(i).class),'double') || strcmp((w(i).class),'single')
            coordsOrDistSize = w(i).size;
            if coordsOrDistSize(2) == 2 && coordsOrDistSize(1) > 2
                coordinates{i} = w(i);
            elseif coordsOrDistSize(1) == coordsOrDistSize(2) && coordsOrDistSize(1) > 2 
                distances{i} = w(i);
            end
        end
    end
    
    coordinates(cellfun(@isempty,coordinates)) = [];
    distances(cellfun(@isempty,distances)) = [];
    
    if isequal(app.tabOtherCoordsOrDist.SelectedObject,app.tabOtherCoordinates) 
        if ~isempty(coordinates)
            names = arrayfun(@(i)coordinates{i}.name,1:size(coordinates,1),'UniformOutput',false);
            app.tabOtherChooseFromWs.String = ['Choose from workspace';names'];
        else
            app.tabOtherChooseFromWs.String = 'Choose from workspace';
        end
        
        app.tabOtherdistType.String = 'Distance Type: EUC';
    elseif isequal(app.tabOtherCoordsOrDist.SelectedObject,app.tabOtherDistance)
        if ~isempty(distances)
            names = arrayfun(@(i)distances{i}.name,1:size(distances,1),'UniformOutput',false);
            app.tabOtherChooseFromWs.String = ['Choose from workspace';names'];
        else
            app.tabOtherChooseFromWs.String = 'Choose from workspace';
        end
        app.tabOtherdistType.String = 'Distance Type: EXPLICIT';
    else
        app.tabOtherChooseFromWs.String = 'Choose from workspace';
    end
    app.tabOtherChooseFromWs.Value = 1;
    app.tabOthernCities.String = 'Number of cities: - - -';

end
