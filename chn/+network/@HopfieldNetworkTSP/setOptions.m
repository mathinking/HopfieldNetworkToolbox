function net = setOptions(net,opts)
    
    net.SimFcn   = opts.SimFcn;
	net.TrainFcn = opts.TrainFcn;
 
    assert(opts.TrainParam.K <= floor(net.TrainParam.N/2), 'HopfieldNetworkTSP:invalid_value', ...
        'Input ''K'' (number of chains) must be less or equal than floor(N/2)');
    net.TrainParam.K = opts.TrainParam.K;
    
    if ~isempty(opts.Cities.Names)
        assert(length(opts.Cities.Names) == net.TrainParam.N, 'HopfieldNetworkTSP:invalid_value', ...
            ['City names must match number of neurons (', num2str(net.TrainParam.N), ').']);
        net.Cities.Names = opts.Cities.Names; 
    else
        net.Cities.Names = net.cityTextGeneration(net.TrainParam.N);
    end

    net.Cities.DistanceType = opts.Cities.DistanceType;
    
    % Errors (not assertions) should have been taken care in
    % tsphopfieldnetOptions
    if isempty(opts.Cities.Coordinates)
        if isempty(opts.Cities.DistanceMatrix) % No Coordinates or Distance Matrix provided
            if ~strcmp(opts.Cities.DistanceType, 'explicit') 
                net.Cities.Coordinates = net.polygonCoords(1, net.TrainParam.N);
                net.Cities.DistanceMatrix = net.computeDistance(net.Cities.Coordinates, net.Cities.DistanceType);    
            else
                error('HopfieldNetworkTSP:invalid_value','''DistanceMatrix'' is empty but ''DistanceType'' property is set to ''explicit''.');
            end
        else % Coordinates not provided, Distance Matrix provided
            assert(size(opts.Cities.DistanceMatrix,1) == net.TrainParam.N, 'HopfieldNetworkTSP:invalid_value', ...
                ['''DistanceMatrix'' must have ', num2str(net.TrainParam.N), ' rows and columns']);            
            if ~strcmp(opts.Cities.DistanceType, 'explicit')
                error('HopfieldNetworkTSP:DistanceMatrixGivenTypeNotExplicit', ...
                    'Set ''DistanceType'' to ''explicit'' if ''DistanceMatrix'' is provided and ''Coordinates'' are empty.');
            else
                net.Cities.Coordinates = opts.Cities.Coordinates; % Empty
                net.Cities.DistanceMatrix = opts.Cities.DistanceMatrix; 
                if ~isfloat(net.Cities.DistanceMatrix)
                    warning('HopfieldNetworkTSP:setOptions:DataConversion', ['''DistanceMatrix'' converted from ', class(net.Cities.DistanceMatrix), ' to double']);
                    net.Cities.DistanceMatrix = double(net.Cities.DistanceMatrix);
                end
            end
        end
    else
        assert(size(opts.Cities.Coordinates,1) == net.TrainParam.N, 'HopfieldNetworkTSP:invalid_value', ...
            ['City ''Coordinates'' matrix must match the number of neurons: ', num2str(net.TrainParam.N)]);

        if isempty(opts.Cities.DistanceMatrix) % Coordinates provided, Distance Matrix not provided
            if ~strcmp(opts.Cities.DistanceType, 'explicit')
                net.Cities.Coordinates = opts.Cities.Coordinates;
                if ~isfloat(net.Cities.Coordinates)
                    warning('HopfieldNetworkTSP:setOptions:DataConversion', ['''Coordinates'' converted from ', class(net.Cities.Coordinates), ' to double']);
                    net.Cities.Coordinates = double(net.Cities.Coordinates);
                end
                net.Cities.DistanceMatrix = net.computeDistance(net.Cities.Coordinates, net.Cities.DistanceType);    
            else
                error('HopfieldNetworkTSP:invalid_value','''DistanceMatrix'' is empty but ''DistanceType'' property is set to ''explicit''.');
            end
        else % Coordinates provided, Distance Matrix provided
            assert(size(opts.Cities.DistanceMatrix,1) == net.TrainParam.N, 'HopfieldNetworkTSP:invalid_value', ...
                ['''DistanceMatrix'' must have ', num2str(net.TrainParam.N), ' rows and columns']);            
            if ~strcmp(opts.Cities.DistanceType, 'explicit')
                net.Cities.Coordinates = opts.Cities.Coordinates;
                if ~isfloat(net.Cities.Coordinates)
                    warning('HopfieldNetworkTSP:setOptions:DataConversion', ['''Coordinates'' converted from ', class(net.Cities.Coordinates), ' to double']);
                    net.Cities.Coordinates = double(net.Cities.Coordinates);
                end
                net.Cities.DistanceMatrix = opts.Cities.DistanceMatrix;
                if ~isfloat(net.Cities.DistanceMatrix)
                    warning('HopfieldNetworkTSP:setOptions:DataConversion', ['''DistanceMatrix'' converted from ', class(net.Cities.DistanceMatrix), ' to double']);
                    net.Cities.DistanceMatrix = double(net.Cities.DistanceMatrix);
                end
            else
                error('HopfieldNetworkTSP:nonEmptyCoordinates', ...
                    '''Coordinates'' must be empty if ''DistanceType'' is set to ''explicit''.');                 
            end            
        end        
    end
    
    if ~isempty(opts.Cities.Subtours) 
        net.verifyIfValidSubtours(opts.Cities.Subtours, opts.Cities.SubtoursPositions, net.Cities.Names);
    end
    net.Cities.Subtours = opts.Cities.Subtours;
    net.Cities.SubtoursPositions = opts.Cities.SubtoursPositions;
    
    if ~isempty(opts.Cities.Tau)
        assert(opts.Cities.Tau <= net.TrainParam.N - 1, 'HopfieldNetworkTSP:invalid_value', 'Tau must be less or equal than N - 1');
        net.Cities.Tau = opts.Cities.Tau;
    end
    if isempty(opts.Cities.Tau) && strcmp(opts.SimFcn,'divide-conquer')
        net.Cities.Tau = max(3,round(net.TrainParam.N/10));
    end    
    if isempty(opts.Cities.Tau) && ~strcmp(opts.SimFcn,'divide-conquer')
        net.Cities.Tau = net.TrainParam.N - 1;
    end

    net.Cities = orderfields(net.Cities);
 
end
