function net = init(net)
    
	net.Results.Time = zeros(1,net.Setting.MaxIter);
    net.Results.Energy = nan(1,net.Setting.MaxIter); 
    net.Results.Energy(1) = 0;
    net.Results.ExitFlag = [];
    net.Results.CompTime = [];
    net.Results.ItersReached = [];

end
