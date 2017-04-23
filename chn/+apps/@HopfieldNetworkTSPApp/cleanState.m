function cleanState(app)

    app.net = [];

    if isequal(app.tabgroup.SelectedTab,app.tabOther)  
        cb_tabOtherCoordsDist(app)
    end

    defaultSettings(app);

    app.settings_u0Edit.Enable = 'off';
    app.settings_dtEdit.Enable = 'off';
    app.settings_maxIterEdit.Enable = 'off';
    app.settings_eEdit.Enable = 'off';
    app.settings_R_ITEREdit.Enable = 'off';
    app.settings_qEdit.Enable = 'off';
    app.cities_tauEdit.Enable = 'off';

    app.elapsedTime.String = '';

    delete(app.plot.Children)
    delete(app.energyplot.Children)
    app.plot.XLim = [0,1];
    app.plot.YLim = [0,1];
    app.plot.XTick = 0:0.2:1;
    app.plot.YTick = 0:0.1:1;
    app.plot.YScale = 'linear';
    app.plot.Title.String = 'Problem Coordinates';

    app.energyplot.XLim = [0,1];
    app.energyplot.YLim = [0,1];
    app.energyplot.XTick = 0:0.2:1;
    app.energyplot.YTick = 0:0.2:1;
    app.energyplot.YScale = 'linear';

end
