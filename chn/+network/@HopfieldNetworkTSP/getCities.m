function returningValue = getCities(net,field)   

    if nargin == 1;
        returningValue = net.cities;
    else
        if isfield(net.cities,field)
            returningValue = net.cities.(field);
        else
            error('hopfieldNetwork:unvalidSetting',['Unvalid setting: ', field]);
        end
    end
    
end
