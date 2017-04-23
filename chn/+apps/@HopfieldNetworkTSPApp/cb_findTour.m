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

if app.simFcnMenu.Value == 1
    msgbox('Please choose an algorithm for the Hopfield Network simulation', 'Error', 'error');
	return;
else
    algorithms = app.simFcnMenu.String;
    setSimFcn(app.net,algorithms{app.simFcnMenu.Value});
    if strcmp(algorithms{app.simFcnMenu.Value},'divide-conquer')
        setCities(app.net,'tau',max(3,round(getTrainParam(app.net,'N')/10))); %#TODO Bring to APP
    end
end

u0 = str2double(app.settings_u0Edit.String);
maxIter = str2double(app.settings_maxIterEdit.String);
R_ITER = str2double(app.settings_R_ITEREdit.String);
dt = str2double(app.settings_dtEdit.String);
e = str2double(app.settings_eEdit.String);
q = str2double(app.settings_qEdit.String);

try
    setSetting(app.net,'u0',u0)
    setSetting(app.net,'maxIter',maxIter)
    setSetting(app.net,'R_ITER',R_ITER)
    setSetting(app.net,'dt',dt)
    setSetting(app.net,'e',e)
    setSetting(app.net,'q',q)
catch me
	msgbox(me.message)
    return;
end

if isequal(app.hwResources.SelectedObject,app.hwResourcesCPU)
	setSetting(app.net,'hwResources','CPU')
else
	setSetting(app.net,'hwResources','GPU')
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

if getResults(app.net, 'exitFlag') ~= 1
    msgbox('You may need to increase the number of iterations','Warning','warn');
else
    plot(app.net,'total',[],app.plot);
    energyplot(app.net,app.energyplot);    
end
