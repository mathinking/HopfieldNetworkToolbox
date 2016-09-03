function disp(net)
    if length(net) > 1
        fprintf([num2str(size(net,1)),'x',num2str(size(net,2)), ...
            ' hopfield array \n\n']);
        return;
    end

    tab1 = 13;
    tab2 = 21;
    strSection = @(rightTab,section) ...
        ['<strong>',repmat(' ',1,rightTab-length(section)),...
        section,': ','</strong>'];
    fprintf([strSection(tab1,'Neural Network'),' ','\n\n']);
    fprintf([strSection(tab1,'name'),net.name,'\n\n']);
    fprintf([strSection(tab1,'trainFcn'),net.trainFcn,'\n']);
	fprintf('\n');
    fprintf([strSection(tab1,'simFcn'),net.simFcn,'\n']);
    fprintf('\n');
    fprintf([strSection(tab1,'trainParam'),'\n']);
    dispStruct(net.trainParam);
    fprintf('\n');
    fprintf([strSection(tab1,'originalParameters'),'\n']);
    dispStruct(net.originalParameters);
    fprintf('\n');
	fprintf([strSection(tab1,'problemParameters'),'\n']);
    dispStruct(net.problemParameters);
    fprintf('\n');
    fprintf([strSection(tab1,'setting'),'\n']);
    dispStruct(net.setting);
    if ~isempty(net.results.x)
        fprintf('\n');
        fprintf([strSection(tab1,'results'),'\n']);
        dispStruct(net.results);
    end
    fprintf('\n');

    function dispStruct(thisStruct)
        structfields = fields(thisStruct);
        for i = 1:length(structfields)
            thisfield = thisStruct.(structfields{i});
            if isscalar(thisfield) && isnumeric(thisfield)
                field2str = num2str(thisfield);
            elseif ischar(thisfield)
                field2str = thisfield;
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
    end
end
