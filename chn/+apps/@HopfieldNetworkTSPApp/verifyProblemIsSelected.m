function problemIsSelected = verifyProblemIsSelected(app)

    problemIsSelected = true;
    if isequal(app.tabgroup.SelectedTab,app.tabTSPLIB) && app.tabTSPLIBmenu.Value == 1
        problemIsSelected = false;
    elseif isequal(app.tabgroup.SelectedTab,app.tabPolygon) && isnan(str2double(app.tabPolygonnCitiesEdit.String))
        problemIsSelected = false;
    elseif isequal(app.tabgroup.SelectedTab,app.tabOther) && strcmp(app.tabOtherChooseFromWs.String{app.tabOtherChooseFromWs.Value},'Choose from workspace')
        problemIsSelected = false;
    end

end