function returningValue = getSetting(net,field)
    if nargin == 1;
        returningValue = net.setting;
    else
        if isfield(net.setting,field)
            returningValue = net.setting.(field);
        else
            error('hopfieldNetwork:unvalidSetting',['Unvalid setting: ', field]);
        end
    end
end
