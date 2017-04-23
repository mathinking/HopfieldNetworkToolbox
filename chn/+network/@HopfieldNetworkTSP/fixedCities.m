function [net,V,U] = fixedCities(net)

    myCities = net.cities.names;
    dMax = max(max(net.cities.d));
    V = zeros(net.trainParam.N);
    citiesInTours = cellfun(@(C)strsplit(C,'-'),net.cities.fixedCities,'UniformOutput',false);
    citiesInTours = [citiesInTours{:}];
    freeCities = setdiff(myCities,citiesInTours); %freeCities
    [~,freeCitiesNumber] = intersect(myCities,freeCities);

    for f = 1:length(net.cities.fixedCities)
        % Modify distance matrix
        thisTour = strsplit(net.cities.fixedCities{f},'-');

        % Fat Points.
        % Distances for cities belonging to a fat point are
        % modified to 0
        [~,~,tourOrder] = intersect(thisTour,myCities,'stable');
        net.cities.d(tourOrder,tourOrder) = 0;

        % Avoid free cities join an already 2-cities connected city
        % Free cities (cities not connected to any other city)
        % cannot connect to cities that already have 2 connections.
        middleTour = thisTour(2:end-1);
        if ~isempty(middleTour)     
            [~,~,middleTourNumber] = intersect(middleTour,myCities,'stable');   
            net.cities.d(middleTourNumber,freeCitiesNumber) = 10*dMax;
            net.cities.d(freeCitiesNumber,middleTourNumber) = 10*dMax;   
        end
        % All other distances remain the same, as connections are
        % feasible.

        % Starting Point V
        % Starting point V is automatically modified to take into
        % account all connections that have been fixed.
        for c = 1:length(tourOrder)
            V(tourOrder(c),tsphopfieldnet.modulo(net.cities.startFixedCitiesIn(f) + c-1,net.trainParam.N)) = 1;
        end
    end
    % Neurons whose cities and positions have not been fixed are
    % randomly created.
    randRows = sum(V,2) == 0;
    randCols = sum(V,1) == 0;
    V(randRows,randCols) = rand(nnz(randRows),nnz(randCols)) * 1e-5;
    U = net.setting.invTransferFcn(V);

end