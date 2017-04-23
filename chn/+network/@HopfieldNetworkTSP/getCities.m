function returningValue = getCities(net,field)   

    if nargin == 1
        returningValue = net.Cities;
    else
        if isfield(net.Cities,field)
            returningValue = net.Cities.(field);
        else
            error('hopfieldNetwork:unvalidSetting',['Unvalid setting: ', field]);
        end
    end
    
end
