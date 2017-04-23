function net = init(net)
    
    net = init@network.HopfieldNetwork(net);
    
    net.Results.ValidPath = false;
    net.Results.TourLength = [];
    net.Results.VisitOrder = [];
    
    net.Results = orderfields(net.Results);
    
end
