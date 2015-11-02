function net = init(net)
    
    net = init@hopfieldnetwork(net);
    
    net.results.validSolution = false;
    net.results.x = [];
    
    net.results = orderfields(net.results);
    
end
