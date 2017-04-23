function verifyIfValidSubtours(subtours, startingPositions, myCities)
    
    if any(cellfun(@isempty,subtours))
        error('HopfieldNetworkTSP:invalidsubtours','Please make sure that there are no empty Subtours.');
    end    
    if any(startingPositions > length(myCities))
        error('HopfieldNetworkTSP:invalidsubtoursPositions','Please make sure you have provided valid starting positions (less than the number of cities) for each Subtour.');
    end
    assert(length(startingPositions) == length(subtours), 'HopfieldNetworkTSP:unevensubtours', ...
        'The number of Subtours starting positions does not match the number of Subtours.');
        
    N = length(myCities);

    fixedPositions = zeros(length(subtours),N);

    for s = 1:length(subtours)
        thisSubtourCities = strsplit(subtours{s},'-');
        try
            cityPos = cellfun(@(c) find(strcmp(myCities,c)),thisSubtourCities);
        catch 
            error('HopfieldNetworkTSP:SubtourCityNotPresentInNames','Verify that the cities forming the Subtour are present in the ''Names'' property.');
        end

        [n, bin] = histc(cityPos, unique(cityPos));
        multiple = find(n > 1);
        repeatedCitiesPos = find(ismember(bin, multiple));

        if ~isempty(repeatedCitiesPos)
            if ~(length(repeatedCitiesPos) == 2 && (repeatedCitiesPos(1) == 1 && repeatedCitiesPos(end) == N+1))
                repeatedCities = unique(thisSubtourCities(repeatedCitiesPos));
                if length(repeatedCities) == 1
                    error('HopfieldNetworkTSP:SubtourCityRepeatedSameSubtour',['City ',repeatedCities{1},' is repeated in subtour ',subtours{s}]);
                else
                    error('HopfieldNetworkTSP:SubtourCitiesRepeatedSameSubtour',['Cities ', strjoin(repeatedCities,', '),' are repeated in subtour ',subtours{s}]);
                end
            end
        end

        startPos = startingPositions(s);

        fixedPositions(s,cityPos) = network.HopfieldNetworkTSP.modulo(startPos:startPos+length(thisSubtourCities)-1,N);

    end

    k = find(any(fixedPositions));

    for p = k
        thisCityPositions = fixedPositions(fixedPositions(:,p) ~= 0,p);

        if length(unique(thisCityPositions)) > 1
            positionsStr = sprintf('%d and ',sort(thisCityPositions));
            positionsStr(end-4:end) = [];
            error('HopfieldNetworkTSP:verifyIfValidSubtours:repeatedCity',['City ',myCities{p},' is repeated in different subtours in positions: ', positionsStr,'.'])
        end
    end
    
    % Repeated Subtour Positions
    fixedPositions = sum(fixedPositions);
    fixedPositions = fixedPositions(fixedPositions > 0);
    assert(length(unique(fixedPositions)) == length(fixedPositions), ...
        'HopfieldNetworkTSP:verifyIfValidSubtours:repeatedSubtourPosition', ...
        'There are repeated positions in provided Subtours.');
end
