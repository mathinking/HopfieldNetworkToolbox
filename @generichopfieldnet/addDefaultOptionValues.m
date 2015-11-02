function net = addDefaultOptionValues(net, options)

%     if isempty(options.cities.names)
%         net.cities.names = tsphopfieldnet.cityTextGeneration(net.trainParam.N);
%     end
% 
%     if isempty(options.cities.type)
%         net.cities.type = 'EUC';
%     end
% 
%     if isempty(options.cities.coords) && isempty(options.cities.d)
%         net.cities.coords = tsphopfieldnet.polygonCoords(1, net.trainParam.N);
%     end
%     
%     if isempty(options.cities.d)
%     	net.cities.d = [];
%     end
% 
%     if isempty(options.cities.fixedCities)
%         net.cities.fixedCities = {''};
%         net.cities.startFixedCitiesIn = NaN;
%     end
    
    if isempty(options.trainFcn)
        net.trainFcn = 'traingty';
    end
    
    if isempty(options.simFcn)
        net.simFcn = 'talavan-yanez';
    end
        
end
