function cb_parameterCEdit(app,~)

C = str2double(app.parameterCEdit.String);
if isempty(C)
    msgbox('You must provide a value for C.', 'Error', 'error');
    return;
end
if isnan(C)
    msgbox('Provide a numeric value for C', 'Error', 'error');
    app.parameterCEdit.String = '';
    return;
end
if C <= 0
    msgbox('Provide a positive value for C', 'Error', 'error');
    return;
end

app.elapsedTime.String = '';

% Only update if net already exists
if isequal(app.tabgroup.SelectedTab,app.tabTSPLIB)
    if app.tabTSPLIBmenu.Value > 1 
        updateTspHopfieldNet(app,C);
    end
elseif isequal(app.tabgroup.SelectedTab,app.tabPolygon) 
    if ~isempty(app.tabPolygonnCitiesEdit.String)
        updateTspHopfieldNet(app,C);
    end
elseif isequal(app.tabgroup.SelectedTab,app.tabOther)
    if app.tabOtherChooseFromWs.Value > 1
        updateTspHopfieldNet(app,C);
    end
end
