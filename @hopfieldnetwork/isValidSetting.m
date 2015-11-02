function isValidSetting(property,value)
    
    invalidDataTypeID = 'hopfieldnetwork:invalid_datatype';
    
    switch property
        case{'u0','maxIter','e','q','R_ITER','dt','viewConvergenceSpeed'}
            assert(isa(value,'double'), invalidDataTypeID, [property, ' must be double.']);
        case{'hwResources','transferFcn'}
            assert(isa(value,'char'), invalidDataTypeID, [property, ' must be char.']);
        case{'showCommandLine','loggingV','viewConvergence'}
            assert(isa(value,'logical'), invalidDataTypeID, [property, ' must be logical.']);
    end
    
    invalidValueID = 'hopfieldnetwork:invalid_value';
    
    switch property
        case 'u0'
            assert(value > 0, invalidValueID, '''u0'' must be greater than 0.');
            
        case 'hwResources'
            assert(strcmp(value,'CPU') || strcmp(value,'GPU'), invalidValueID, '''hwResources'' must be ''CPU'' or ''GPU'''); 
 
        case 'maxIter'
            assert(value > 0, invalidValueID, 'Maximum number of iterations must be greater than 0.');
            
        case 'e'
            assert(value > 6 & value < 22, invalidValueID, 'Exponent ''e'' must be between 6 and 22.');
            
        case 'q'
            assert(value > 0 & value < 1, invalidValueID, '''q'' must be between 0 and 1.');
            
        case 'R_ITER'
            assert(value >= 0, invalidValueID, 'R_ITER must be greater or equal to 0.');
            
        case 'dt'
            assert(value > 0, invalidValueID, '''dt'' must be greater than 0.');
                       
        case 'transferFcn'
            assert(strcmp(value,'satlin') || strcmp(value,'tanh'), invalidValueID, 'Transfer function must be ''satlin'' or ''tanh''');

        case 'invTransferFcn'
            assert(strcmp(value,'invsatlin') || strcmp(value,'atanh'), invalidValueID, 'Inverse transfer function must be ''invsatlin'' or ''atanh'''); 
            
        case 'viewConvergenceSpeed'
            assert(value > 0 & value < 1 , invalidValueID, '''viewConvergenceSpeed'' must be between 0 and 1.');
    end
end
