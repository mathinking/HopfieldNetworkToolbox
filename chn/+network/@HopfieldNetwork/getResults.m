function returningValue = getResults(net,field)
    if nargin == 1
        returningValue = net.Results;
    else
        if isfield(net.Results,field)
            returningValue = net.Results.(field);
        else
            error('hopfieldNetwork:unvalidSetting',['Unvalid setting: ', field]);
        end
    end
end
