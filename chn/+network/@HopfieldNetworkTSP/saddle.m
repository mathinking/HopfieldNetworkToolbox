function S = saddle(net, method)

	simFcn = net.simFcn;
    if ~strcmp(simFcn,'talavan-yanez')
        net.simFcn = 'talavan-yanez';
    end    

    if nargin < 2 
        method = 'numeric';
    end

    if strcmp(method, 'numeric')
        V = ones(net.trainParam.N,net.trainParam.N-net.trainParam.K)/(net.trainParam.N-net.trainParam.K);
        U = net.setting.invTransferFcn(V);
        if ~isfield(net.trainParam, 'Np')
            net = train(net);
        end
        S = sim(net,V,U);
        reinit(net);
        
    elseif strcmp(method, 'numeric2') %#TODO Make this the default method
        if ~isfield(net.trainParam, 'Np')
            net = train(net);
        end        
        vxi_col = (((net.trainParam.N + 1)*net.trainParam.C + 3) * ...
            ones(net.trainParam.N) + ((net.trainParam.N-2)*(net.trainParam.C + 3) ...
            -(net.trainParam.N-1)*net.trainParam.rho) * eye(net.trainParam.N) + ...
            2*net.trainParam.D*net.cities.d) \ ones(net.trainParam.N,1) * net.trainParam.C * net.trainParam.Np;
        S = repmat(vxi_col,1,net.trainParam.N);
              
    elseif strcmp(method, 'analytic')

        A = net.trainParam.A;
        B = net.trainParam.B;
        C = net.trainParam.C;
        D = net.trainParam.D;
        N = net.trainParam.N;
        Np = net.trainParam.Np;
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
        net.simFcn = simFcn;
    end
end
