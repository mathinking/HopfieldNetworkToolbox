function net = setOptions(net,options)

    if ~isempty(options.cities.coords)
        assert(size(options.cities.coords,1) == net.trainParam.N && size(options.cities.coords,2) == 2, 'tsphopfieldnet:invalid_value', ...
            ['City coordinates matrix ''coords'' must have ', num2str(net.trainParam.N), ' rows and 2 columns']);
    end
    
	net.cities.coords = options.cities.coords; 

    if ~isempty(options.cities.names)
        assert(length(options.cities.names) == net.trainParam.N, 'tsphopfieldnet:invalid_value', ...
            ['City names must match number of neurons (', num2str(net.trainParam.N), ').']);
        assert(all(cellfun(@ischar,options.cities.names)), 'tsphopfieldnet:invalid_value', ...
            'City names must be a cell array of strings.');
    end

    net.cities.names = options.cities.names;
    net.cities.type = options.cities.type;
    net.cities.fixedCities = options.cities.fixedCities;
    net.cities.startFixedCitiesIn = options.cities.startFixedCitiesIn;
    net.cities.tau = options.cities.tau;
    net.trainParam.K = options.trainParam.K;
    % Verify here fixedCities.
    
    if ~isempty(options.cities.d)
        assert(size(options.cities.d,1) == net.trainParam.N && size(options.cities.d,2) == net.trainParam.N, 'tsphopfieldnet:invalid_value', ...
            ['Matrix distance ''d'' must have ', num2str(net.trainParam.N), ' rows and columns']);
        net.cities.d = options.cities.d;

    elseif ~strcmp(options.cities.type,'EXPLICIT') && ~isempty(options.cities.type)
        net.cities.d = tsphopfieldnet.computeDistance(net.cities.coords,net.cities.type);

    end
    
    net.simFcn = options.simFcn;
	net.trainFcn = options.trainFcn;
 
end
