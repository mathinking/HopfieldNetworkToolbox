function returningValue = getTrainParam(net,field)
    if nargin == 1
        returningValue = net.TrainParam;
    else
        if isfield(net.TrainParam,field)
            returningValue = net.TrainParam.(field);
        else
            error('hopfieldNetwork:unvalidSetting',['Unvalid setting: ', field]);
        end
    end
end
