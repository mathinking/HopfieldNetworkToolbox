function net = addDefaultOptionValues(net, options)

    if isempty(options.cities.names)
        net.cities.names = tsphopfieldnet.cityTextGeneration(net.trainParam.N);
    end

    if isempty(options.cities.type)
        net.cities.type = 'EUC';
    end

    if isempty(options.cities.coords) && isempty(options.cities.d)
        net.cities.coords = tsphopfieldnet.polygonCoords(1, net.trainParam.N);
    end
    
    if isempty(options.cities.d)
    	net.cities.d = [];
    end

    if isempty(options.cities.fixedCities)
        net.cities.fixedCities = {''};
        net.cities.startFixedCitiesIn = NaN;
    end
    
    if isempty(options.trainFcn)
        net.trainFcn = 'trainty';
    end
    
    if isempty(options.simFcn)
        net.simFcn = 'talavan-yanez';
    end

    if isempty(options.trainParam.K)
        net.trainParam.K = 0;
    end
    
    if isempty(options.cities.tau) && strcmp(options.simFcn,'divide-conquer')
        net.cities.tau = max(3,round(net.trainParam.N/10));
    end
    
    if isempty(options.cities.tau) && (isempty(options.simFcn) || strcmp(options.simFcn,'euler') || strcmp(options.simFcn,'talavan-yanez'))
        net.cities.tau = net.trainParam.N - 1;
    end    
    
end
