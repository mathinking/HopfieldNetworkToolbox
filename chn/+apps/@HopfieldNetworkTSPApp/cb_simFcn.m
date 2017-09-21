function cb_simFcn(app,~)

    problemIsSelected = verifyProblemIsSelected(app);
    if ~problemIsSelected
        msgbox('Please select first the problem to solve', 'Error', 'error');
        app.simFcnMenu.Value = 1;
        return;
    end
    if app.schemeMenu.Value == 1
        msgbox('Please select first the Hopfield Network scheme to solve the problem', 'Error', 'error');
        app.simFcnMenu.Value = 1;
        return;
    end
    
    algorithms = app.simFcnMenu.String;

    chosenAlgorithm = algorithms{app.simFcnMenu.Value};

    app.elapsedTime.String = '';
       
    if strcmp(chosenAlgorithm,'euler') || strcmp(chosenAlgorithm,'runge-kutta')
        app.settings_u0Edit.Enable = 'on';
        app.settings_dtEdit.Enable = 'on';
        app.settings_maxIterEdit.Enable = 'on';
        app.settings_eEdit.Enable = 'on';
        app.settings_R_ITEREdit.Enable = 'on';
        app.settings_qEdit.Enable = 'on';
        if ~isempty(app.net)
            setSimFcn(app.net,chosenAlgorithm)
        end

    elseif strcmp(chosenAlgorithm,'talavan-yanez')
        app.settings_u0Edit.Enable = 'on';
        app.settings_dtEdit.Enable = 'off';
        app.settings_maxIterEdit.Enable = 'on';
        app.settings_eEdit.Enable = 'on';
        app.settings_R_ITEREdit.Enable = 'on';
        app.settings_qEdit.Enable = 'on';
        if ~isempty(app.net)
            setSimFcn(app.net,'talavan-yanez')
        end
    else
        app.settings_u0Edit.Enable = 'off';
        app.settings_dtEdit.Enable = 'off';
        app.settings_maxIterEdit.Enable = 'off';
        app.settings_eEdit.Enable = 'off';
        app.settings_R_ITEREdit.Enable = 'off';
        app.settings_qEdit.Enable = 'off';
    end
   
end
