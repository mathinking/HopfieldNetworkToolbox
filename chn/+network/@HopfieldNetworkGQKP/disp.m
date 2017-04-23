function disp(net)
    if length(net) > 1
        fprintf([num2str(size(net,1)),'x',num2str(size(net,2)), ...
            ' hopfieldnet array \n\n']);
        return;
    end
    tab = 27;
    strSection = @(rightTab,section) ...
        [repmat(' ',1,rightTab-length(section)),section,': '];
    %fprintf([strSection(tab,['  ',net.Name,' with properties']),' ','\n\n']);
    
    fprintf([repmat(' ',1,tab-length(['  ',net.Name,' with properties'])),...
        '  ','<a href = "matlab:helpPopup ',net.Name, '">', net.Name,...
        '</a>',' with properties',': ','\n\n']);    
    printsection('TrainFcn', net);
	fprintf('\n');
    printsection('SimFcn', net);
    fprintf('\n');
    printsection('TrainParam');
    dispStruct(net.TrainParam);
    fprintf('\n');
    printsection('OriginalParameters');
    dispStruct(net.OriginalParameters);
    fprintf('\n');
	printsection('ProblemParameters');
    dispStruct(net.ProblemParameters);
    fprintf('\n');
    printsection('Setting'); 
    dispStruct(net.Setting);
    if ~isempty(net.Results.x)
        fprintf('\n');
        printsection('Results');
        dispStruct(net.Results);
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

            fprintf([strSection(tab,structfields{i}),...
                field2str,'\n']);
        end       
    end
    function printsection(field, net)
        if nargin == 1
            fprintf([strSection(tab,field),'\n']);
        else
            fprintf([strSection(tab,field),net.(field),'\n']);
        end
        fprintf([repmat(' ',1,tab-length(field)),repmat('=',1,length(field)),'\n']);
    end
end
