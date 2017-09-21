function cb_tabPolygonNcities(app)
   
    N = str2double(app.tabPolygonnCitiesEdit.String);
    if isnan(N) || N < 2
        msgbox('You must provide a value for number of cities greater than 2', 'Error', 'error');
        return;
    end
    createTspHopfieldNet(app,N)
    
end

