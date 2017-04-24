function setCities(net,property,value)

    if isfield(net.Cities,property)
        
        switch property
            case 'DistanceType'
                % Cannot change from explicit to any other DistanceType
                % without using tsphopfieldnetOptions
                assert(~strcmp(net.Cities.DistanceType,'explicit'),'tsphopfieldnet:invalid_value',...
                    ['Trying to modify an ''explicit'' distance matrix. You may provide a new ''DistanceType''',...
                    ' together with new ''Coordinates'' using the tsphopfieldnetOptions class.']);
                if ~strcmp(value,'explicit') % Recompute distance matrix. 
%                     net.Cities.DistanceMatrix = tsphopfieldnet.computeDistance(net.Cities.Coordinates, value);
                    options = tsphopfieldnetOptions('DistanceType',value,'DistanceMatrix',net.Cities.DistanceMatrix);
                    net.Cities.DistanceType = options.Cities.DistanceType;
                else
                    net.Cities.DistanceType = value;
                end
            case 'DistanceMatrix'
                options = tsphopfieldnetOptions('DistanceMatrix',value,'DistanceType','explicit');
                assert(size(value,1) == net.TrainParam.N, 'tsphopfieldnet:invalid_value', ...
                    ['''DistanceMatrix'' must have ', num2str(net.TrainParam.N), ' rows and columns']);
                net.Cities.DistanceMatrix = options.Cities.DistanceMatrix;

            case 'Subtours'
                options = tsphopfieldnetOptions('Subtours',value,'SubtoursPositions',net.Cities.SubtoursPositions,...
                    'Names',net.Cities.Names);
                net.Cities.Subtours = options.Cities.Subtours;
                
            case 'SubtoursPositions'
                options = tsphopfieldnetOptions('SubtoursPositions',value,'Subtours',net.Cities.Subtours,...
                    'Names',net.Cities.Names);
                net.Cities.SubtoursPositions = options.Cities.SubtoursPositions;
            
            case 'Names'
                options = tsphopfieldnetOptions('Names',value,'Subtours',net.Cities.Subtours,...
                    'SubtoursPositions',net.Cities.SubtoursPositions);
                assert(length(options.Cities.Names) == net.TrainParam.N, 'tsphopfieldnet:invalid_value', ...
                    ['City names must match number of neurons (', num2str(net.TrainParam.N), ').']);                
                net.Cities.Names = options.Cities.Names;
            
            case 'Coordinates'
                % Recompute distance matrix
                options = tsphopfieldnetOptions(property,value,'DistanceType',net.Cities.DistanceType);
                assert(size(options.Cities.(property),1) == net.TrainParam.N, 'tsphopfieldnet:invalid_value', ...
                    ['City ''Coordinates'' matrix must match the number of neurons: ', num2str(net.TrainParam.N)]);
                
            case 'Tau'
                options = tsphopfieldnetOptions(property,value,'SimFcn',net.SimFcn);
                if ~isempty(options.Cities.Tau)
                    assert(options.Cities.Tau <= net.TrainParam.N - 1, 'tsphopfieldnet:invalid_value', 'Tau must be less or equal than N - 1');
                    net.Cities.Tau = options.Cities.Tau;
                end
                if isempty(options.Cities.Tau) && strcmp(options.SimFcn,'divide-conquer')
                    net.Cities.Tau = max(3,round(net.TrainParam.N/10));
                end
                if isempty(options.Cities.Tau) && ~strcmp(options.SimFcn,'divide-conquer')
                    net.Cities.Tau = net.TrainParam.N - 1;
                end
        end
       
    else
        error('hopfieldNetwork:unvalidSetting',['There is no Cities property named ', property]);
    end

end
