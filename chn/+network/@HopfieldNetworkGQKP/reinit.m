function net = reinit(net)

    if ~isempty(net.Results.x) % Re-init results
        net = init(net);
    end
    
end
