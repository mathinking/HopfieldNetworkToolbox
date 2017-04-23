function S = saddle(net, method)

	simFcn = net.SimFcn;
    if ~strcmp(simFcn,'talavan-yanez')
        net.SimFcn = 'talavan-yanez';
    end    

    if nargin < 2 
        method = 'numeric';
    end

    if strcmp(method, 'numeric')
        V = ones(net.TrainParam.N,net.TrainParam.N-net.TrainParam.K)/(net.TrainParam.N-net.TrainParam.K);
        U = net.Setting.InvTransferFcn(V);
        if ~isfield(net.TrainParam, 'Np')
            net = train(net);
        end
        S = sim(net,V,U);
        reinit(net);
        
    elseif strcmp(method, 'numeric2') %#TODO Make this the default method
        if ~isfield(net.TrainParam, 'Np')
            net = train(net);
        end        
        vxi_col = (((net.TrainParam.N + 1)*net.TrainParam.C + 3) * ...
            ones(net.TrainParam.N) + ((net.TrainParam.N-2)*(net.TrainParam.C + 3) ...
            -(net.TrainParam.N-1)*net.TrainParam.rho) * eye(net.TrainParam.N) + ...
            2*net.TrainParam.D*net.cities.d) \ ones(net.TrainParam.N,1) * net.TrainParam.C * net.TrainParam.Np;
        S = repmat(vxi_col,1,net.TrainParam.N);
              
    elseif strcmp(method, 'analytic')

        A = net.TrainParam.A;
        B = net.TrainParam.B;
        C = net.TrainParam.C;
        D = net.TrainParam.D;
        N = net.TrainParam.N;
        Np = net.TrainParam.Np;
        d = net.cities.d;

        T = repmat({zeros(N)},N,N);
        % Some additional matrices needed for creating T
        M2 = ones(N);
        I = eye(N);

        % Adding sigma Matrix
        %d = d + sigma_fun_2(N,1);
        for i = 1:N
            T{i,i} = A*(M2-I);
        end

        pos = find(~I)';
        for i = pos
            T{i} = T{i} + B*I;
        end

        for i = 1:N^2
            T{i} = T{i} + C*M2;
        end

        [JJ,II] = meshgrid(1:N,1:N);
        aux = mod(JJ,N) == mod(II+1,N) | mod(JJ,N) == mod(II-1,N);
        for x = 1:N
            for i = 1:N
                T{x,i} = T{x,i} + D*d(x,i)*aux;
            end
        end

        M = cell2mat(T);
        b = Np*C*ones(N^2,1);

        S = reshape(M\b,N,N)';

    else
        error('tsphopfieldnet:saddle', 'Unvalid saddle method calculation.');
    end
    
    if ~strcmp(simFcn,'talavan-yanez')
        net.SimFcn = simFcn;
    end
end
