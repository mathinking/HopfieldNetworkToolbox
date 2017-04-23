function net = init(net)
    
	net.results.time = zeros(1,net.setting.maxIter);
    net.results.energy = nan(1,net.setting.maxIter); 
    net.results.energy(1) = 0;
    net.results.exitFlag = [];
    net.results.compTime = [];
    net.results.itersReached = [];

end