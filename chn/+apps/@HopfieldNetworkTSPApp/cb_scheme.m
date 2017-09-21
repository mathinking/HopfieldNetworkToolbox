function cb_scheme(app,~)

    problemIsSelected = verifyProblemIsSelected(app);
    if ~problemIsSelected
        msgbox('Please select first the problem to solve', 'Error', 'error');
        app.schemeMenu.Value = 1;
        return;
    end
    schemes = app.schemeMenu.String;

    chosenScheme = schemes{app.schemeMenu.Value};

    if strncmp(chosenScheme,'divide-conquer',14)
        app.cities_tauEdit.Enable = 'on';
        tau = max(3,round(getTrainParam(app.net,'N')/10));
        setCities(app.net,'Tau',tau); 
        app.cities_tauEdit.String = num2str(tau);
        
        train(app.net);
        trainParams = getTrainParam(app.net);

        app.parameterAEdit.String = num2str(trainParams.A);
        app.parameterBEdit.String = num2str(trainParams.B);
        app.parameterDEdit.String = num2str(1/trainParams.dUaux);
        
        app.parameterNEdit.String = num2str(trainParams.N);
        app.parameterNpEdit.String = num2str(trainParams.Np);
        app.cities_tauEdit.String = num2str(getCities(app.net,'Tau'));        
    else
        app.cities_tauEdit.Enable = 'off';
    end
    
    if ~isempty(app.net) && ~strcmp(chosenScheme,'Scheme')
        setScheme(app.net, chosenScheme)
    end

end
