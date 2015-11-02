function isValidSetting(property,value)
    
    invalidDataTypeID = 'tsphopfieldnet:invalid_datatype';
    
    switch property
        case {'d','coords','startFixedCitiesIn','tau','K'}
            assert(isa(value,'double'), invalidDataTypeID, [property, ' must be double.']);
            
        case {'names','fixedCities'}
            assert(iscell(value), invalidDataTypeID, [property, ' must be cell array of strings.']);
            
        case 'type'
            assert(isa(value,'char'), invalidDataTypeID, [property, ' must be char.']);
            
        case {'trainFcn', 'simFcn'}
            assert(isa(value,'char'), invalidDataTypeID, [property, ' must be char.']);
    end
    
    invalidValueID = 'tsphopfieldnet:invalid_value';
    
    switch property
        case 'd'
            assert(size(value,1) == size(value,2), invalidValueID, [property, ' must be a square matrix']);
            
        case 'coords'
            assert(size(value,2) == 2 || isempty(value), invalidValueID, [property, ' must be a data point matrix with two columns']);
            
        case 'trainFcn'
            assert(strcmp('trainty',value), invalidValueID, 'Training function (parametrization) must be ''trainty''');
            
        case 'simFcn'
            assert(strcmp('euler',value) || strcmp('talavan-yanez',value) || strcmp('divide-conquer',value), invalidValueID, 'Simulation function must be ''euler'' or ''talavan-yanez'' or ''divide-conquer''');
            
    end
end
