function setTrainParam(net,property,value)

    if isfield(net.TrainParam,property)
        net.TrainParam.(property) = value;
    else
        error(['There is no training parameter named ', property]);
    end
    
end
