function net = reinit(net)

    if ~isempty(net.Results.TourLength) % Re-init results
        net = init(net);
    end
    
end
