function cb_findTour(app,~)

if isequal(app.tabgroup.SelectedTab,app.tabTSPLIB)
    if app.tabTSPLIBmenu.Value == 1
        msgbox('Please choose a TSPLIB problem', 'Error', 'error');
        return;
    end
elseif isequal(app.tabgroup.SelectedTab,app.tabPolygon)
    N = str2double(app.tabPolygonnCitiesEdit.String);
    if isnan(N) 
        msgbox('You must provide a value for number of cities', 'Error', 'error');
        return;
    end
elseif isequal(app.tabgroup.SelectedTab,app.tabOther)
    if app.tabOtherChooseFromWs.Value == 1
        msgbox('You must create a variable from the workspace and select it from the list below', 'Error', 'error');
        return;        
    end
end

if app.schemeMenu.Value == 1
    msgbox('Please choose a scheme for the Hopfield Network simulation', 'Error', 'error');
	return;
end

u0 = str2double(app.settings_u0Edit.String);
maxIter = str2double(app.settings_maxIterEdit.String);
R_ITER = str2double(app.settings_R_ITEREdit.String);
dt = str2double(app.settings_dtEdit.String);
e = str2double(app.settings_eEdit.String);
q = str2double(app.settings_qEdit.String);
tau = str2double(app.cities_tauEdit.String);

try
    setSetting(app.net,'U0',u0)
    setSetting(app.net,'MaxIter',maxIter)
    setSetting(app.net,'R_Iter',R_ITER)
    setSetting(app.net,'Dt',dt)
    setSetting(app.net,'E',e)
    setSetting(app.net,'Q',q)
    setCities(app.net,'Tau',tau)
catch me
	msgbox(me.message)
    return;
end

if isequal(app.ExecutionEnvironment.SelectedObject,app.ExecutionEnvironmentCPU)
	setSetting(app.net,'ExecutionEnvironment','CPU')
else
	setSetting(app.net,'ExecutionEnvironment','GPU')
end

if isequal(app.seed.SelectedObject,app.seedShuffle)
	rng('shuffle');
else
    seed = str2double(app.seedFixedEdit.String);
    if isnan(seed) || seed <= 0
        msgbox('Provide a positive value for the seed')
    end
    rng(seed);
end

start(app.elapsedTimer);

train(app.net);
sim(app.net);
stop(app.elapsedTimer);

if getResults(app.net, 'ExitFlag') ~= 1
    msgbox('You may need to increase the number of iterations','Warning','warn');
else
    plot(app.net,'total',[],app.plot);
    energyplot(app.net,app.energyplot);    
end
