function net = addDefaultOptionValues(net, options)

    if isempty(options.setting.u0)
        net.setting.u0 = 0.3;
    end

    if isempty(options.setting.transferFcn)
        net.setting.transferFcn = @(u) 0.5 * (1 + tanh(u ./ net.setting.u0));
    end   
    
    if isempty(options.setting.invTransferFcn)
        net.setting.invTransferFcn = @(v) (net.setting.u0 / 2) * log(v./(1-v));
    end

    if isempty(options.setting.showCommandLine)
        net.setting.showCommandLine = true;
    end

    if isempty(options.setting.hwResources)
        net.setting.hwResources = 'CPU';
    end

    if isempty(options.setting.maxIter)
        net.setting.maxIter = 2000;
    end

    if isempty(options.setting.e)
        net.setting.e = 13;
    end

    if isempty(options.setting.q)
        net.setting.q = 0.8;
    end

    if isempty(options.setting.R_ITER)
        net.setting.R_ITER = 20;
    end

    if isempty(options.setting.dt)
        net.setting.dt = 0.01; % For Euler's method
    end

    if isempty(options.setting.loggingV)
        net.setting.loggingV = false; 
    end

    if isempty(options.setting.viewConvergence)
        net.setting.viewConvergence = false;
        net.setting.viewConvergenceSpeed = NaN; 
    else
        if isempty(options.setting.viewConvergenceSpeed)
            net.setting.viewConvergenceSpeed = 0.3; 
        end
    end

end
