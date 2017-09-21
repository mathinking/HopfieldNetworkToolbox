function cb_closeApp(app,~)

    selection = questdlg('Close Hopfield Net TSP solver App?',...
        'Close App?',...
        'Yes','No','Yes'); 
    switch selection
        case 'Yes'
            warning('on','tsphopfieldnet:NotSimulated');
            warning('on','MATLAB:MKDIR:DirectoryExists');
            delete(app.figure);
        case 'No'
            return;
    end
end