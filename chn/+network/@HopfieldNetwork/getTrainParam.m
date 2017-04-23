function returningValue = getTrainParam(net,field)
    if nargin == 1;
        returningValue = net.trainParam;
    else
        if isfield(net.trainParam,field)
            returningValue = net.trainParam.(field);
        else
            error('hopfieldNetwork:unvalidSetting',['Unvalid setting: ', field]);
        end
    end
end