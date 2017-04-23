function setTrainParam(net,property,value)

    if isfield(net.TrainParam,property)
        if strcmp(property, 'N')
            error('tsphopfieldnet:setTrainParam:NonModifiable', ...
                'Parameter ''N'' cannot be modified. Rebuild the network with the ''tsphopfieldnet'' constructor.');
        elseif strcmp(property, 'K')
            options = tsphopfieldnetOptions('K', value);   
            assert(options.TrainParam.K <= floor(net.TrainParam.N/2), 'tsphopfieldnet:invalid_value', ...
                'Input ''K'' (number of chains) must be less or equal than floor(N/2)');
            net.SimFcn = options.SimFcn;
        end
            
        net.TrainParam.(property) = value;
    else
        error('tsphopfieldnet:setTrainParam:invalid_parameter',['There is no training parameter named ', property]);
    end
    
end
