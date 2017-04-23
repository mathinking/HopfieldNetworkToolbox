function options = createOptions(varargin)  

    options.cities.coords = []; 
    options.cities.names = {}; 
    options.cities.d = [];
    options.cities.type = '';
    options.cities.fixedCities = [];
    options.cities.startFixedCitiesIn = [];
    options.cities.tau = [];
    
    options.simFcn = '';
    options.trainFcn = '';
    
    options.trainParam.K = [];
    
    numargs = nargin;

    i = 1; extra = 1;
    while i < numargs
        field = varargin{i};
        value = varargin{i+1};
        [options,hopfieldnetoptions{extra}] = setfield(options,field,value); %#ok<AGROW>
            i = i + 2;
            if ~isempty(hopfieldnetoptions{extra})
                extra = extra + 1;
            end
    end
    
    if exist('hopfieldnetoptions','var') && ~isempty(hopfieldnetoptions{1})
        hopfieldnetoptions = [hopfieldnetoptions{:}];
        options2 = hopfieldnetwork.createOptions(hopfieldnetoptions{:});
        options.setting = options2.setting;
    end
    
    if numargs == 0 || ~isfield(options,'setting')
        options2 = hopfieldnetwork.createOptions();
        options.setting = options2.setting;        
    end
    
end

function [options,hopfieldnetoptions] = setfield(options,field,value)

    hopfieldnetoptions = [];
    switch field
        case {'d','coords','names','type','fixedCities','startFixedCitiesIn','tau'}
            tsphopfieldnet.isValidSetting(field,value);
            options.cities.(field) = value;
            
        case {'trainFcn'}
            tsphopfieldnet.isValidSetting(field,value);
            options.trainFcn = value;
            
        case {'simFcn'}
            tsphopfieldnet.isValidSetting(field,value);
            options.simFcn = value;
            
        case {'K'}
            tsphopfieldnet.isValidSetting(field,value);
            options.trainParam.K = value;
            
        otherwise % Is an option defined in hopfieldnetwork class
            hopfieldnetoptions = {field,value};
    end

end
