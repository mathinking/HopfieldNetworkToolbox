function S = saddle(net, method)

	simFcn = net.SimFcn;
    if ~strcmp(simFcn,'talavan-yanez')
        net.SimFcn = 'talavan-yanez';
    end    

    if nargin < 2 
        method = 'analytic';
    end

    if strcmp(method, 'analytic')
        
        if ~isfield(net.TrainParam, 'Np')
            net = train(net);
        end        
        
        One = ones(net.TrainParam.N);
        I = eye(net.TrainParam.N);

        JK = I;
        for i = 1:net.TrainParam.K
            JK(2*i-1:2*i,2*i-1:2*i) = ones(2) - eye(2);
        end
        DK = JK * net.Cities.DistanceMatrix;
        DK(logical(I)) = 0;

        IK = (net.TrainParam.K > 0) * JK + I;
        vxi_col = ((net.TrainParam.B + (net.TrainParam.N - net.TrainParam.K)*net.TrainParam.C) * One + ...
            (net.TrainParam.N - net.TrainParam.K - 1) * net.TrainParam.A * IK - net.TrainParam.B*I + ...
            net.TrainParam.D * (DK + DK')) \ ones(net.TrainParam.N,1) * net.TrainParam.C * net.TrainParam.Np;

        S = repmat(vxi_col,1,net.TrainParam.N-net.TrainParam.K);
        
    elseif strcmp(method, 'numeric') 
        
        V = ones(net.TrainParam.N,net.TrainParam.N-net.TrainParam.K)/(net.TrainParam.N-net.TrainParam.K);
        U = net.Setting.InvTransferFcn(V);
        if ~isfield(net.TrainParam, 'Np')
            net = train(net);
        end
        S = sim(net,V,U,true);
        reinit(net);
              
    else
        error('tsphopfieldnet:saddle', 'Unvalid saddle method calculation.');
    end
    
    if ~strcmp(simFcn,'talavan-yanez')
        net.SimFcn = simFcn;
    end
end
