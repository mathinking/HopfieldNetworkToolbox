function disp(net)
    if length(net) > 1
        fprintf([num2str(size(net,1)),'x',num2str(size(net,2)), ...
            ' hopfield array \n\n']);
        return;
    end

    tab = 27;
    strSection = @(rightTab,section) ...
        [repmat(' ',1,rightTab-length(section)),section,': '];
%     fprintf([strSection(tab,['  ',net.Name,' with properties']),' ','\n\n']);
    
    fprintf([repmat(' ',1,tab-length(['  ',net.Name,' with properties'])),...
        '  ','<a href = "matlab:helpPopup ',net.Name, '">', net.Name,...
        '</a>',' with properties',': ','\n\n']);
    printsection('TrainFcn',net);
	fprintf('\n');
    printsection('SimFcn',net)
    fprintf('\n');
    printsection('TrainParam')
    dispStruct(net.TrainParam);
    fprintf('\n');
    printsection('Cities')
    dispStructCities(net.Cities);
    fprintf('\n');
    printsection('Setting')
    dispStructSetting(net.Setting);
    if ~isempty(net.Results.TourLength)
        fprintf('\n');
        printsection('Results')
        dispStruct(net.Results);
    end
    fprintf('\n');

    function dispStructCities(thisStructCities)
        structCitiesFields = fields(thisStructCities);
        for i = 1:length(structCitiesFields)
            thisfieldName = structCitiesFields{i};
            thisfieldValues = thisStructCities.(thisfieldName);
            
            if strcmp(thisfieldName, 'Subtours') || strcmp(thisfieldName, 'Names')
                if ~isempty(thisfieldValues)
                    field2str = strjoin(thisfieldValues, ', ');
                    maxField2str = 20;
                    commas = strfind(field2str,',');

                    if any(commas > maxField2str)
                        firstPositionDelete = find(commas > maxField2str,1,'first');
                        numMissingNames = length(commas(firstPositionDelete:end));
                        field2str(commas(firstPositionDelete):end) = [];
                        if numMissingNames == 1
                            if strcmp(thisfieldName, 'Subtours')
                                endingstr = ' more subtour';                            
                            else
                                endingstr = ' more city';
                            end
                        else
                            if strcmp(thisfieldName, 'Subtours')
                                endingstr = ' more subtours';
                            else
                                endingstr = ' more cities';
                            end
                        end
                        field2str = [field2str, ' and ', num2str(numMissingNames),endingstr]; %#ok<AGROW>
                    end 

                    fprintf([strSection(tab,structCitiesFields{i}),field2str,'\n']);      
                else
                    printfield(structCitiesFields, thisfieldValues, i, strSection, tab);
                end
                
            elseif strcmp(thisfieldName, 'SubtoursPositions')
                printfield(structCitiesFields, thisfieldValues, i, strSection, tab);

            elseif strcmp(thisfieldName, 'Tau')
                printfield(structCitiesFields, thisfieldValues, i, strSection, tab);
            else
                printfield(structCitiesFields, thisfieldValues, i, strSection, tab);
            end                
        end
    end
    function dispStruct(thisStruct)
        structfields = fields(thisStruct);
        for i = 1:length(structfields)
            thisfield = thisStruct.(structfields{i});
            printfield(structfields, thisfield, i, strSection, tab);
        end       
    end
    function dispStructSetting(thisStructSetting)
        structfields = fields(thisStructSetting);
        for i = 1:length(structfields)
            thisfieldName = structfields{i};
            thisfieldValues = thisStructSetting.(thisfieldName);
            if strcmp(thisfieldName, 'CheckpointPath')
                thisfieldValues = regexprep(thisfieldValues,'\\','\\\\');
            end
            printfield(structfields, thisfieldValues, i, strSection, tab);
        end   
    end
    function printsection(field, net)
        if nargin == 1
            fprintf([strSection(tab,field),'\n']);
        else
            fprintf([strSection(tab,field),net.(field),'\n']);
        end
        fprintf([repmat(' ',1,tab-length(field)),...
            repmat('=',1,length(field)),'\n']);
    end
end

function printfield(structfields, thisfield, i, strSection, tab2)
    
    if isnumeric(thisfield) && isempty(thisfield)
        field2str = '[]';
    elseif isscalar(thisfield) && isnumeric(thisfield)
        field2str = num2str(thisfield);
    elseif ischar(thisfield)
        field2str = thisfield;
        if isempty(thisfield)
            field2str = '''''';
        end
    elseif isnumeric(thisfield)
        field2str = ['[',num2str(size(thisfield,1)),...
            'x',num2str(size(thisfield,2)),' ', ...
            num2str(class(thisfield)),']'];
    elseif isa(thisfield,'function_handle')
        field2str = func2str(thisfield);
    elseif islogical(thisfield)
        if thisfield
            field2str = 'true';
        else
            field2str = 'false';
        end
    end    

    fprintf([strSection(tab2,structfields{i}),...
        field2str,'\n']);
end
