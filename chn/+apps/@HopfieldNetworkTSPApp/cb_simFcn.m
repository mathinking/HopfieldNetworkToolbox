function cb_simFcn(app,~)

    algorithms = app.simFcnMenu.String;

    chosenAlgorithm = algorithms(app.simFcnMenu.Value);

    app.elapsedTime.String = '';

    if strcmp(chosenAlgorithm,'euler')
        app.settings_u0Edit.Enable = 'on';
        app.settings_dtEdit.Enable = 'on';
        app.settings_maxIterEdit.Enable = 'on';
        app.settings_eEdit.Enable = 'on';
        app.settings_R_ITEREdit.Enable = 'on';
        app.settings_qEdit.Enable = 'on';
        app.cities_tauEdit.Enable = 'off';
        if ~isempty(app.net)
            setSimFcn(app.net,'euler')
        end

    elseif strcmp(chosenAlgorithm,'talavan-yanez')
        app.settings_u0Edit.Enable = 'on';
        app.settings_dtEdit.Enable = 'off';
        app.settings_maxIterEdit.Enable = 'on';
        app.settings_eEdit.Enable = 'on';
        app.settings_R_ITEREdit.Enable = 'on';
        app.settings_qEdit.Enable = 'on';
        app.cities_tauEdit.Enable = 'off';
        if ~isempty(app.net)
            setSimFcn(app.net,'talavan-yanez')
        end    

    elseif strcmp(chosenAlgorithm,'divide-conquer')
        app.settings_u0Edit.Enable = 'on';
        app.settings_dtEdit.Enable = 'off';
        app.settings_maxIterEdit.Enable = 'on';
        app.settings_eEdit.Enable = 'on';
        app.settings_R_ITEREdit.Enable = 'on';
        app.settings_qEdit.Enable = 'on';
        app.cities_tauEdit.Enable = 'on';
        if ~isempty(app.net)
            setSimFcn(app.net,'divide-conquer')
        end 

    else
        app.settings_u0Edit.Enable = 'off';
        app.settings_dtEdit.Enable = 'off';
        app.settings_maxIterEdit.Enable = 'off';
        app.settings_eEdit.Enable = 'off';
        app.settings_R_ITEREdit.Enable = 'off';
        app.settings_qEdit.Enable = 'off';
        app.cities_tauEdit.Enable = 'off';
    end

end
