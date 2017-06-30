function S = saddle(net, method)

    method = 'analytic'; % No other method currently supported for HopfieldNetworkGQKP
    
    if ~isfield(net.TrainParam, 'T')
        net = train(net);
    end        
    S = -net.TrainParam.T\net.TrainParam.ib;    
    
end