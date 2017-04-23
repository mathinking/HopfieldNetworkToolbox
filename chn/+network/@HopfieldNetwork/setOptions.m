function net = setOptions(net, options)
   
    net.Setting.U0 = options.Setting.U0;
    net.Setting.Verbose = options.Setting.Verbose;
    net.Setting.VerboseFrequency = options.Setting.VerboseFrequency;
    net.Setting.ExecutionEnvironment = options.Setting.ExecutionEnvironment;
    net.Setting.MaxIter = options.Setting.MaxIter;
    net.Setting.E = options.Setting.E;
    net.Setting.Q = options.Setting.Q;
    net.Setting.R_Iter = options.Setting.R_Iter;
    net.Setting.Dt = options.Setting.Dt;
    net.Setting.CheckpointPath = options.Setting.CheckpointPath;
    net.Setting.SimulationPlot = options.Setting.SimulationPlot;
    net.Setting.SimulationPlotPauseTime = options.Setting.SimulationPlotPauseTime;
    
    % Transfer function must be either satlin or tanh
    if strcmp(options.Setting.TransferFcn,'satlin')
        net.Setting.TransferFcn = @(u) net.satlin(u, net.Setting.U0);
        net.Setting.InvTransferFcn = @(v) net.invsatlin(v, net.Setting.U0);
    else
        net.Setting.TransferFcn = @(u) 0.5 * (1 + tanh(u ./ net.Setting.U0));
        net.Setting.InvTransferFcn = @(v) (net.Setting.U0 / 2) * log(v./(1-v));
    end
    
end
