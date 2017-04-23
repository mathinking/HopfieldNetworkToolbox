function net = init(net)
    
    net = init@network.HopfieldNetwork(net);
    
    net.Results.ValidSolution = false;
    net.Results.x = [];
    
    net.Results = orderfields(net.Results);
    
end
