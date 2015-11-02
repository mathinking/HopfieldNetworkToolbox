function isValidSetting(property,value)
    
    invalidDataTypeID = 'generichopfieldnet:invalid_datatype';
    
    switch property
%         case {'d','coords','startFixedCitiesIn'}
%             assert(isa(value,'double'), invalidDataTypeID, [property, ' must be double.']);
            
%         case {'names','fixedCities'}
%             assert(iscell(value), invalidDataTypeID, [property, ' must be cell array of strings.']);
            
%         case 'type'
%             assert(isa(value,'char'), invalidDataTypeID, [property, ' must be char.']);
            
        case {'trainFcn', 'simFcn'}
            assert(isa(value,'char'), invalidDataTypeID, [property, ' must be char.']);
    end
    
    invalidValueID = 'generichopfieldnet:invalid_value';
    
    switch property
%         case 'd'
%             assert(size(value,1) == size(value,2), invalidValueID, [property, ' must be a square matrix']);
%             
%         case 'coords'
%             assert(size(value,2) == 2 || isempty(value), invalidValueID, [property, ' must be a data point matrix with two columns']);
            
        case 'trainFcn'
            assert(strcmp('traingty',value), invalidValueID, 'Training function (parametrization) must be ''traingty''');
            
        case 'simFcn'
            assert(strcmp('euler',value) || strcmp('talavan-yanez',value), invalidValueID, 'Simulation function must be ''euler'' or ''talavan-yanez''');
            
    end
end
