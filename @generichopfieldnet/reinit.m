function net = reinit(net)

    if ~isempty(net.results.x) % Re-init results
        net = init(net);
    end
    
end
