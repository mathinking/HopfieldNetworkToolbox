function options = createOptions(varargin)

    options.setting.R_ITER = [];
    options.setting.dt = []; % For Euler's method
    options.setting.e = [];
    options.setting.hwResources = '';
    options.setting.maxIter = [];
    options.setting.q = [];
    options.setting.showCommandLine = [];
    options.setting.transferFcn = '';
    options.setting.invTransferFcn = '';    
    options.setting.u0 = [];
    options.setting.loggingV = [];
    options.setting.viewConvergence = [];
    options.setting.viewConvergenceSpeed = [];

    numberargs = nargin;

    i = 1;
    while i < numberargs
        field = varargin{i};
        value = varargin{i+1};
        options = setfield(options,field,value);       
        i = i + 2;
    end
end

%--------------------------------------------------------------------------
function options = setfield(options,field,value)

	hopfieldnetwork.isValidSetting(field,value);

    switch field
        case {'u0','hwResources','maxIter','e','q','R_ITER','dt',...
                'showCommandLine','loggingV','viewConvergence',...
                'viewConvergenceSpeed'}
            options.setting.(field) = value;

        case {'transferFcn','invTransferFcn'}
            if strcmp(value,'satlin') || strcmp(value,'invsatlin') 
                options.setting.transferFcn = 'satlin';
                options.setting.invTransferFcn = 'invsatlin';
            else
                options.setting.transferFcn = 'tanh';
                options.setting.invTransferFcn = 'atahn';
            end

        otherwise
            error('hopfieldNetwork:unvalidSetting',['Unvalid setting: ', field]);    
    end
end
