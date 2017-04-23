function returningValue = getSetting(net,field)
    if nargin == 1
        returningValue = net.Setting;
    else
        if isfield(net.Setting,field)
            returningValue = net.Setting.(field);
        else
            error('hopfieldNetwork:unvalidSetting',['Unvalid setting: ', field]);
        end
    end
end
