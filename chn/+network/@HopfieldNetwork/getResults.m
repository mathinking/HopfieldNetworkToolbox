function returningValue = getResults(net,field)
    if nargin == 1;
        returningValue = net.results;
    else
        if isfield(net.results,field)
            returningValue = net.results.(field);
        else
            error('hopfieldNetwork:unvalidSetting',['Unvalid setting: ', field]);
        end
    end
end