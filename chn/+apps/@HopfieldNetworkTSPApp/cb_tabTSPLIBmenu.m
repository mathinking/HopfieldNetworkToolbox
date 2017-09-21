function cb_tabTSPLIBmenu(app)

problems = app.tabTSPLIBmenu.String;
problemSelected = app.tabTSPLIBmenu.Value;

app.schemeMenu.Value = 1;
app.simFcnMenu.Value = 1;
app.cities_tauEdit.Enable = 'off';

% Useful for problems that take longer to compute tsplib
% and training
app.tabTSPLIBnCitiesEdit.String = '- - - ';
app.tabTSPLIBdistTypeEdit.String = '- - -';
app.parameterAEdit.String = '';
app.parameterBEdit.String = '';
app.parameterDEdit.String = '';
app.parameterNEdit.String = '';
app.parameterNpEdit.String = '';
app.elapsedTime.String = '';

drawnow;
if problemSelected > 1
    problem = tsplib(problems(problemSelected));
    app.tabTSPLIBnCitiesEdit.String = num2str(problem.NumberOfCities);
    app.tabTSPLIBdistTypeEdit.String = problem.DistanceType;
    createTspHopfieldNet(app,problem)
end
