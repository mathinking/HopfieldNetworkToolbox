function net = init(net)
    
    net = init@hopfieldnetwork(net);
    
    net.results.validPath = false;
    net.results.tourLength = [];
    net.results.visitOrder = [];
    
    net.results = orderfields(net.results);
    
end
