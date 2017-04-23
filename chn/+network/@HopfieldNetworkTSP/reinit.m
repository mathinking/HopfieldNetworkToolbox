function net = reinit(net)

    if ~isempty(net.results.tourLength) % Re-init results
        net = init(net);
    end
    
end