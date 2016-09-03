function isValidSetting(property,value)
    
    invalidDataTypeID = 'generichopfieldnet:invalid_datatype';
    
    switch property            
        case {'trainFcn', 'simFcn'}
            assert(isa(value,'char'), invalidDataTypeID, [property, ' must be char.']);
    end
    
    invalidValueID = 'generichopfieldnet:invalid_value';
    
    switch property
            
        case 'trainFcn'
            assert(strcmp('traingty',value), invalidValueID, 'Training function (parametrization) must be ''traingty''');
            
        case 'simFcn'
            assert(strcmp('euler',value) || strcmp('talavan-yanez',value), invalidValueID, 'Simulation function must be ''euler'' or ''talavan-yanez''');
            
    end
end
