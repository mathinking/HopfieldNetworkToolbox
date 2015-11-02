function net = verifyIfValidSubtours(net)

    subtours = net.cities.fixedCities;
    startingPositions = net.cities.startFixedCitiesIn;
    myCities = net.cities.names;
    
    if length(net.cities.fixedCities) ~= length(net.cities.startFixedCitiesIn)
        error('tsphopfieldnet:unevensubtours','Please provide a starting position for all subtours');
    end
    if any(cellfun(@isempty,subtours))
        error('tsphopfieldnet:invalidsubtour','Please make sure that there are no empty subtours');
    end    

    N = length(myCities);

    fixedPositions = zeros(length(subtours),N);

    for s = 1:length(subtours)
        thisSubtourCities = strsplit(subtours{s},'-');
        try
            cityPos = cellfun(@(c) find(strcmp(myCities,c)),thisSubtourCities);
        catch 
            error('Check that the cities forming the subtour are present in the problem');
        end

        [n, bin] = histc(cityPos, unique(cityPos));
        multiple = find(n > 1);
        repeatedCitiesPos = find(ismember(bin, multiple));

        if ~isempty(repeatedCitiesPos)
            if ~(length(repeatedCitiesPos) == 2 && (repeatedCitiesPos(1) == 1 && repeatedCitiesPos(end) == N+1))
                repeatedCities = unique(thisSubtourCities(repeatedCitiesPos));
                if length(repeatedCities) == 1
                    error(['City ',repeatedCities{1},' is repeated in subtour ',subtours{s}]);
                else
                    error(['Cities ',repeatedCities{:},' are repeated in subtour ',subtours{s}]);
                end
            end
        end

        startPos = startingPositions(s);

        fixedPositions(s,cityPos) = tsphopfieldnet.modulo(startPos:startPos+length(thisSubtourCities)-1,N);

    end

    k = find(any(fixedPositions));

    for p = k
        thisCityPositions = fixedPositions(fixedPositions(:,p) ~= 0,p);

        if length(unique(thisCityPositions)) > 1
            positionsStr = sprintf('%d and ',sort(thisCityPositions));
            positionsStr(end-4:end) = [];
            error(['City ',myCities{p},' is repeated in different subtours in positions: ', positionsStr,'.'])
        end
    end

end
