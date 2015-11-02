function net = setOptions(net, options)
   
    net.setting.R_ITER = options.setting.R_ITER;
    net.setting.dt = options.setting.dt; % For Euler's method
    net.setting.e = options.setting.e;
    net.setting.hwResources = options.setting.hwResources;
    net.setting.maxIter = options.setting.maxIter;
    net.setting.showCommandLine = options.setting.showCommandLine;
    net.setting.q = options.setting.q;
    net.setting.u0 = options.setting.u0;
    net.setting.loggingV = options.setting.loggingV;
    net.setting.viewConvergence = options.setting.viewConvergence;
    net.setting.viewConvergenceSpeed = options.setting.viewConvergenceSpeed;

    % Transfer function must be either satlin or tanh
    if strcmp(options.setting.transferFcn,'satlin')
        net.setting.transferFcn = @(u) net.satlin(u, net.setting.u0);
        net.setting.invTransferFcn = @(v) net.invsatlin(v, net.setting.u0);
    else
        net.setting.transferFcn = @(u) 0.5 * (1 + tanh(u ./ options.setting.u0));
        net.setting.invTransferFcn = @(v) (net.setting.u0 / 2) * log(v./(1-v));
    end
    
%     net.setting.R_ITER = [];
%     net.setting.dt = []; % For Euler's method
%     net.setting.e = [];
%     net.setting.hwResources = '';
%     net.setting.maxIter = [];
%     net.setting.q = [];
%     net.setting.showCommandLine = [];
%     net.setting.u0 = [];
%     net.setting.loggingV = [];
%     net.setting.viewConvergence = [];
%     net.setting.viewConvergenceSpeed = [];
%     net.setting.transferFcn = '';
%     net.setting.invTransferFcn = '';    
%     
%     setSetting(net, 'R_ITER', options.setting.R_ITER)
%     setSetting(net, 'dt', options.setting.dt)
%     setSetting(net, 'e', options.setting.e)
%     setSetting(net, 'hwResources', options.setting.hwResources)
%     setSetting(net, 'maxIter', options.setting.maxIter)
%     setSetting(net, 'showCommandLine', options.setting.showCommandLine)
%     setSetting(net, 'q', options.setting.q)
%     setSetting(net, 'u0', options.setting.u0)
%     setSetting(net, 'loggingV', options.setting.loggingV)
%     setSetting(net, 'viewConvergence', options.setting.viewConvergence)
%     setSetting(net, 'viewConvergenceSpeed', options.setting.viewConvergenceSpeed)
%     setSetting(net, 'transferFcn', options.setting.transferFcn)
%     setSetting(net, 'invTransferFcn', options.setting.invTransferFcn)    
end
