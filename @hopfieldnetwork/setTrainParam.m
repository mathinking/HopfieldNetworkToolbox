function setTrainParam(net,property,value)

    if isfield(net.trainParam,property)
        net.trainParam.(property) = value;
    else
        error(['There is no training parameter named ', property]);
    end
    
end
