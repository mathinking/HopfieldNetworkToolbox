function V = sim(net,V,U)

	[ST,~] = dbstack('-completenames');
    if ~any(strcmp({ST.name},'saddle'))
        timeID = tic;
    end
    
    net = reinit(net);
    
    if ~isfield(net.trainParam,'Np')
        error('tsphopfieldnet:NotTrained', 'Training has not taken place yet. Use train(net).');
    end
    
    if ~isempty(net.cities.fixedCities{1})
        net = verifyIfValidSubtours(net);
        aux_d = net.cities.d;
        [net,V,U] = fixedCities(net);
    end
    
    if strcmp(net.simFcn,'euler')
        if nargin < 2
            [net,V,~,iter] = simEuler(net);
        else
            [net,V,~,iter] = simEuler(net,V,U);
        end
        
    elseif strcmp(net.simFcn,'talavan-yanez')
        if nargin < 2 && isempty(net.cities.fixedCities{1} )
            [net,V,~,iter] = simTalavanYanez(net);
        else
            [net,V,~,iter] = simTalavanYanez(net,V,U);
        end
        
    elseif strcmp(net.simFcn, 'divide-conquer')
        if nargin < 2
            [net,V,~,iter] = simDivideConquer(net);
        else
            [net,V,~,iter] = simDivideConquer(net,V,U);
        end
        
    elseif strcmp(net.simFcn,'talavan-yanez-n')
        if nargin < 2
            [net,V,~,iter] = simTalavanYanezVarN(net);
        else
            [net,V,~,iter] = simTalavanYanezVarN(net,V,U);
        end
        
    else
        error('tsphopfieldnet:UnknownsimFcn', 'Unknown training algorithm');
    end
    
    if ~isempty(net.cities.fixedCities{1})
        net.cities.d = aux_d;
    end   
    if ~any(strcmp({ST.name},'saddle'))
        net = computeTour(net,V,iter);
        net.results.compTime = toc(timeID);

        % Removing unused energy and time elements.
        net.results.energy = net.results.energy(~isnan(net.results.energy));
        indexRemove = find(diff(net.results.time),1,'last') + 1;
        if length(net.results.time) > indexRemove
            net.results.time(indexRemove:end) = [];
        end
%         net.results.time = [net.results.time(1),net.results.time(net.results.time(2:end) ~= 0)];
    end    
end

% --- Simulation Algorithms --- %

% Euler Algorithm
function [net,V,U,iter] = simEuler(net, V, U) 
    N = net.trainParam.N;
    if nargin == 1
        U = rand(N)-.5;   % TODO Different from Pedro's algorithm
        V = 0.5 + 1e-7*U; % TODO Different from Pedro's algorithm
    end
    if strcmp(net.setting.hwResources,'GPU')
        U = gpuArray(U);
        V = gpuArray(V);
    end
    % Stopping criteria
    stopC1 = power(10, -1   * net.setting.e);
    stopC2 = power(10, -1.5 * net.setting.e);
    maxDiffV = 1;
    unstable = false;

    ib = net.trainParam.C * net.trainParam.Np;

    %%%
    % Summing V rows, columns and all elements for optimal perform.
    sumVcol = sum(V);
    sumVrow = sum(V,2);
    sumV = sum(sumVcol);

    % Memory allocation
    dU = zeros(N);

    iter = 1;

    dt = net.setting.dt;
    net.results.time(iter) = dt;

    while iter <= net.setting.maxIter && (maxDiffV > stopC1 || ...
            (maxDiffV > stopC2 && unstable))

        unstable = false;
        for i = 1:N
            if i == N, deltaPrev = 1; else deltaPrev = i+1; end
            if i == 1, deltaNext = N; else deltaNext = i-1; end
            for x = 1:N
                TV = weightMatrixTimesV(...
                    net, deltaPrev, deltaNext, V, sumVrow, ...
                    sumVcol, sumV, x, i);

                dU(x,i) = TV + ib;

                if (V(x,i) < stopC1 && dU(x,i) > 0) || ...
                        (V(x,i) > 1 - stopC1 && dU(x,i) < 0)
                    unstable = true;
                end
                net.results.energy(iter+1) = ...
                    net.results.energy(iter) + ...
                    0.5* V(x,i) * dU(x,i) - ...
                    ib * V(x,i);
            end
        end

        maxDiffV = 0;
        for i = 1:N
            for x = 1:N
                prevV = V(x,i);
                U(x,i) = U(x,i) + dU(x,i)*net.setting.dt;
                V(x,i) = net.setting.transferFcn(U(x,i));
                if abs(prevV - V(x,i)) > maxDiffV
                    maxDiffV = abs(prevV - V(x,i));
                end
            end          
        end   
        sumVcol = sum(V);
        sumVrow = sum(V,2);
        sumV = sum(sumVcol);       

        net.results.time(iter+1) = net.results.time(iter) + ...
            net.setting.dt;
        iter = iter + 1;            
    end

end

% Talaván-Yáñez Algorithm
function [net,V,U,iter] = simTalavanYanez(net,V,U)
    N = net.trainParam.N;
    K = net.trainParam.K;
    if nargin == 1
        U = rand(N,N-K)-.5;   % TODO Different from Pedro's algorithm
        V = 0.5 + 1e-7*U;     % TODO Different from Pedro's algorithm
    end
    if strcmp(net.setting.hwResources,'GPU')
        U = gpuArray(U);
        V = gpuArray(V);
    end
    % Stopping criteria
    stopC1 = power(10, -1   * net.setting.e);
    stopC2 = power(10, -1.5 * net.setting.e);
    maxDiffV = 1;
    unstable = false;

    u_e = net.setting.invTransferFcn(stopC1); 

	ib = net.trainParam.C * net.trainParam.Np;
    
    %%%
    % Summing V rows, columns and all elements for optimal performance
    sumVcol = sum(V);   % Si
    sumVrow = sum(V,2); % Sx
    sumV = sum(sumVcol);

    % Memory allocation
    dU = zeros(N,N-K);

    iter = 1;

    while iter <= net.setting.maxIter && (maxDiffV > stopC1 || ...
            (maxDiffV > stopC2 && unstable))

        dt = 10^100; % Initial value for dt

        if net.setting.loggingV || net.setting.viewConvergence
            [ST,~] = dbstack('-completenames');
            if ~any(strcmp({ST.name},'saddle'))
                if net.setting.loggingV
                    tsphopfieldnet.loggingV(iter,V,dU);
                end
                if net.setting.viewConvergence
                    if iter == 1
                        fV = viewConvergence(iter,V,net);
                    else
                        fV = viewConvergence(iter,V,net,fV);
                    end
                end
            end
        end

        % Computation of the weight matrix T
        % $$ T_{xi,yj} = -(A*\delta_{x,y}*(1-\delta_{i,j}) + B*(1-\delta_{x,y})*
        % \delta_{i,j} + C - D*d_{x*y} * (\delta_{j,i+1} + \delta_{j,i-1}) $$
        TV = weightMatrixTimesVvectorized(net, V, sumVrow, sumVcol, sumV);
        
        dU = TV + ib;

        dV = 2./net.setting.u0 .* V .* (1-V) .*dU;

%       interiorV = V > 0 & V < 1;
        interiorV = U > u_e & U < -u_e;
        %(V > stopC1) & (V < 1-stopC1);
%       borderV = ~interiorV;

        % Computation of dt             
        % In interior states, $\Delta t$ should not make the state get outside the 
        % interval [0,1]
        criteria1 = interiorV & dV < 0;
        if any(any(criteria1))
            VdV = -V./dV;
            dt = min(dt,min(VdV(criteria1)));
        end                 

        criteria2 = interiorV & dV > 0;
        if any(any(criteria2))
            antVdV = (1-V)./dV;
            dt = min(dt,min(antVdV(criteria2)));
        end              

        criteria3 = (U <= u_e & dU > 0) | (U >= -u_e & dU < 0); %(V == 0 & dU > 0) | (V == 1 & dU < 0);
        unstable = any(any(criteria3));
        %%%                    
        % For extreme states (0 or 1), it must be checked using the potential if 
        % the state is stable or unstable.
        % And _dt_ is computed so that the state stays in [0,1].        
        if unstable
            dt = min(dt, min((abs(U(criteria3)) - u_e) ./ ...
                abs(dU(criteria3))));
        end  

        %%% 
        % $S_{1} = \sum_{i = 1}^{n} \frac{dv_{i}(t)}{dt} 
        % (\sum_{i = j}^{n} T_{i,j}v_{j}(t) + i_{i}^b) = 
        % \sum_{i = 1}^{n} \frac{dv_{i}(t)}{dt} \frac{du_{i}(t)}{dt}$
        S1 = sum(sum(dV.*dU));

        % Computing S2 for optimal dt
        % $S_{2} = - \sum_{i = 1}^{n} \sum_{j = 1}^{n} \frac{dv_{i}(t)}{dt} T_{i,j} 
        % \frac{dv_{j}(t)}{dt}$

        sumVcol = sum(dV);
        sumVrow = sum(dV,2);
        sumV = sum(sumVcol);

        TdV = weightMatrixTimesVvectorized(net, dV, sumVrow, sumVcol, sumV);

        S2 = -sum(sum(dV .* TdV));

        %%% 
        % The appropriate value of $\Delta t$ depends on the term $S_{2}$. For the 
        % case $S_{2} <= 0$, the energy function will decrease for every $\Delta t
        % > 0$ and the greater value $\Delta t$ must be chosen. Otherwise, if
        % $S_{2} > 0$, then the value $\Delta t = S_{2}/S_{1}$ will produce
        % the largest decrease of the energy function. 

        % dt = min(dt,S1/S2);
        sw_optimal = false;
        if S2 > 0                 
            dt = min(dt,S1./S2);
            sw_optimal = true;
        end
        if iter < net.setting.R_ITER && ~sw_optimal
            dt = net.setting.q * dt; % Integration step reduction
        end

        % State update
        % Update of states once optimal $\Delta t$ has been computed. Potential 
        % variables also have to be updated. Note that the treatment is different
        % for interior and border values. 

        Vprev = V; 

        % vi(t) is border Value. vi(t+1) might stay or leave the
        % border.
%                 borderV = (V < stopC1) | (V > 1-stopC1);
        borderV = (U <= u_e) | (U >= -u_e);
        U(borderV) = U(borderV) + dU(borderV).*dt;
        V(borderV) = net.setting.transferFcn(U(borderV));
        U(borderV & U <= u_e) = u_e; V(borderV & U <= u_e) = 0; % Stays in the border
        U(borderV & U >=-u_e) =-u_e; V(borderV & U >=-u_e) = 1; % Stays in the border

        % vi(t) is interior value. vi(t+1) might stay in the
        % interior or reach a border.
        VupdateInInterior = interiorV & ...
            (Vprev + dV .* dt > stopC1) & ...
            (Vprev + dV .* dt < 1-stopC1);
        V(VupdateInInterior) = Vprev(VupdateInInterior) + ...
            dV(VupdateInInterior) .* dt; % Stays in the interior

        VupdateInBorder0 = interiorV & (Vprev + dV.*dt <= stopC1);
        VupdateInBorder1 = interiorV & (Vprev + dV.*dt >= 1-stopC1);

        U(VupdateInBorder0) =  u_e;
        U(VupdateInBorder1) = -u_e;
        V(VupdateInBorder0) =  0;
        V(VupdateInBorder1) =  1;

        maxDiffV = max(max(abs(Vprev-V)));

        sumVcol = sum(V);
        sumVrow = sum(V,2);
        sumV = sum(sumVcol);

        %%%
        % $E(t + \Delta t) = E(t) - S_{1}\Delta t + \frac{1}{2}S_{2}\Delta t^{2}$
        if strcmp(net.setting.hwResources,'GPU')
            S1 = gather(S1);
            S2 = gather(S2);
            dt = gather(dt);
        end
        net.results.energy(iter+1) = net.results.energy(iter) - ...
            S1.*dt + 0.5.*S2.*dt.^2;
        net.results.time(iter+1) = net.results.time(iter) + dt;
        iter = iter + 1;
    end

    if strcmp(net.setting.hwResources,'GPU') && nargout > 1
        V = gather(V);
        U = gather(U);
        iter = gather(iter);
    end
     
end

% Divide and Conquer algorithm
function [net,V,U,iter] = simDivideConquer(net,V,U)
    
    % Determine weather to plot the two phases of the problem. Consider
    % bringing this to createOptions
    plotPhases = false;
    myPlot.myCitiesColorP1 = [0.8,0.8,0.8];
    myPlot.myCitiesTextColor = [0,0,0];
    myPlot.myInsideColorP1 = [0,0,0.8];
    myPlot.myTitleP1 = 'Divide and Conquer. Phase 1 problem (TSP_1^{\tau})';    
    myPlot.myCitiesColorP2 = [1,0,0];
    myPlot.myInsideColorP2 = [0.8,0.8,0.8];
    myPlot.myCitiesTextColor = [0,0,0];
    myPlot.myTitleP2 = 'Divide and Conquer. Phase 2 problem (TSP_2^{*})';
    myPlot.myTitle = 'Divide and Conquer. TSP_1^{\tau} + TSP_2^{*}';

    % The Phase 1 problem. TSP^tau_1
    % 2 possible methods:
    %    a. using tau which uses closes tau neighbours
    %    b. using p which uses distances <= p * net.trainParam.dU + (1-p) * net.trainParam.dL
    
    if net.trainParam.K == 0
        if nargin < 2
            V = saddle(net) + (rand(net.trainParam.N) - 0.5) * 1e-5;
            U = net.setting.invTransferFcn(V);
        end
        [chains,V] = simDivideConquerPhase1(net,V,U,plotPhases,myPlot);
        
        % Reached max number of iterations        
        if net.results.exitFlag == 0
            mkdir('Phase1_InsuficientIterations');
            save(fullfile(pwd,'Phase1_InsuficientIterations',regexprep(datestr(datetime),{':',' ','-'},'_')))
            while net.results.exitFlag <= 0
                V = V + (rand(net.trainParam.N) - 0.5) * 1e-5;
                U = net.setting.invTransferFcn(V);
                [chains,V] = simDivideConquerPhase1(net,V,U,plotPhases,myPlot);
            end
            
        % Possible saddle point reached
        elseif net.results.exitFlag == -1 
            mkdir('Phase1_PossibleSaddlePoint');
            save(fullfile(pwd,'Phase1_PossibleSaddlePoint',regexprep(datestr(datetime),{':',' ','-'},'_')))            
            while net.results.exitFlag <= 0
                V = V + (rand(net.trainParam.N) - 0.5) * 1e-5;
                U = net.setting.invTransferFcn(V);
                [chains,V] = simDivideConquerPhase1(net,V,U,plotPhases,myPlot);
            end
        end
        % TODO Add correct number of iterations
    else
        % Construir chains
        chains = cell(1,net.trainParam.K);
        for c = 1:net.trainParam.K
            chains{c} = [2*c-1,2*c];
        end
    end
    
    if (~isempty(chains) && length(chains{1}) <= net.trainParam.N && net.results.exitFlag == 1) || isempty(chains) && net.results.exitFlag == 1

        % Part 2. 3 possible methods:
        % a. fixing cities
        % b. using new distance
        % c. connecting subtours with greedy
        [netPhase2,V2] = simDivideConquerPhase2(net,chains,plotPhases,myPlot);
        % Building final V
        V = zeros(net.trainParam.N);
        iNew = 1;
        iOld = 1;

        % Reached max number of iterations        
        if netPhase2.results.exitFlag == 0
            mkdir('Phase2_InsuficientIterations');
            save(fullfile(pwd,'Phase2_InsuficientIterations',regexprep(datestr(datetime),{':',' ','-'},'_')))
            while netPhase2.results.exitFlag <= 0
                V2 = V2 + (rand(netPhase2.trainParam.N,netPhase2.trainParam.N-netPhase2.trainParam.K) - 0.5) * 1e-5;
                U = netPhase2.setting.invTransferFcn(V2);
                V2 = sim(netPhase2, V2, U);
            end
            
        % Possible saddle point reached
        elseif netPhase2.results.exitFlag == -1 
            mkdir('Phase2_PossibleSaddlePoint');
            save(fullfile(pwd,'Phase2_PossibleSaddlePoint',regexprep(datestr(datetime),{':',' ','-'},'_')))            
            % Stuck in a saddle point. Perturbation required.
            while netPhase2.results.exitFlag <= 0
                V2 = V2 + (rand(netPhase2.trainParam.N,netPhase2.trainParam.N-netPhase2.trainParam.K) - 0.5) * 1e-5;
                U = netPhase2.setting.invTransferFcn(V2);
                V2 = sim(netPhase2, V2, U);
            end
        end
        % TODO Add correct number of iterations        
        
        while iNew <= net.trainParam.N
            
            thisCity = find(strcmp(net.cities.names,netPhase2.cities.names(netPhase2.results.visitOrder(iOld))));

            V(thisCity,iNew) = 1;
            iNew = iNew + 1;
            if netPhase2.results.visitOrder(iOld) <= 2*netPhase2.trainParam.K
                for c = 1:length(chains)
                    if any(chains{c} == thisCity)
                        break;
                    end
                end
                if rem(netPhase2.results.visitOrder(iOld),2) % Start of chain
                    for thisChain = 2:length(chains{c})
                        V(chains{c}(thisChain),iNew) = 1;
                        iNew = iNew + 1;
                    end
                    iOld = iOld + 1;  
                else % End of chain
                    for thisChain = length(chains{c})-1:-1:1
                        V(chains{c}(thisChain),iNew) = 1;
                        iNew = iNew + 1;
                    end
                    iOld = iOld + 1;                    
                end
            end
            iOld = iOld + 1;            
        end
        U = net.setting.invTransferFcn(V);
        iter = net.results.itersReached + netPhase2.results.itersReached - 1;

        if plotPhases
            if isempty(net.results.tourLength);
                init(net);
                net = computeTour(net,V,iter); % Needed for plot phase1 + phase2 to output correctly
            else
                [~,net.results.visitOrder] = max(V); % Needed for plot phase1 + phase2 to output correctly
            end
            plot(net,'phase2',chains,[],myPlot.myCitiesColorP2,myPlot.myCitiesTextColor,myPlot.myInsideColorP1,myPlot.myTitle);
        end    

        % Removing unused energy and time elements.
        net.results.energy = net.results.energy(~isnan(net.results.energy));
        net.results.time = [net.results.time(1),net.results.time(net.results.time(2:end) ~= 0)];

        net.results.energy = [net.results.energy, netPhase2.results.energy(2:end)];
        net.results.time = [net.results.time, net.results.time(end) + netPhase2.results.time(2:end)];
    else
        iter = net.results.itersReached;
    end
end

function [chains,V] = simDivideConquerPhase1(net,V,U,plotPhases,myPlot)

    p_or_tau = net.cities.tau;
	% Backing up distances
    aux_d = net.cities.d;

    [net.cities.d, neighbours] = tsphopfieldnet.neighbourDistance(net, p_or_tau);

    % Starting point close to saddle point
    if nargin < 2
        V = saddle(net) + (rand(net.trainParam.N) - 0.5) * 1e-5;
        U = net.setting.invTransferFcn(V);
    end
    
    [net,V,U,iter] = simTalavanYanez(net,V,U);

    % Retrieving original distance from backup.    
    net.cities.d = aux_d;
    
    % Finding chains such that go from city x to city y (from
    % position j to position j+1) using the original distance.
   
    net = computeTour(net,V,iter);
    % Extract chains from S.
	chains = {};
            
    if net.results.validPath
        finalTour = net.results.visitOrder;
        
        % X: City matrix. Connecting cities from j to j+1
        X = zeros(net.trainParam.N);
        for x = 1:size(X,1)-1
            X(finalTour(x),finalTour(x+1)) = 1;
        end
        X(finalTour(end),finalTour(1)) = 1;
        S = neighbours & X;
        
        if any(any(S))
            s = 1;
            allChains = false;
            while ~allChains
                doneChain = false;
                y = find(any(S),1,'first');
                x = find(S(:,y));
                S(x,y) = 0;

                thisChain = [y,x];
                while ~doneChain
                    y = x;
                    x = find(S(:,y));
                    if isempty(x)
                        doneChain = true;
                    else
                        thisChain = [thisChain,x]; 
                    end
                    S(x,y) = 0;
                end
                
                doneLookingForGreaterChain = false;
                k = 1;
                
                while ~doneLookingForGreaterChain && k <= length(chains)
                    if any(chains{k} == thisChain(1)) || any(chains{k} == thisChain(end))
                        [~, ia, ib] = intersect(chains{k},thisChain);
                        if ia ~=1
                            chains{k} = fliplr(chains{k}); %#ok<*AGROW>
                        end
                        
                        if ib == 1
                        	thisChain = fliplr(thisChain);
                        end
                        chains{k} = [thisChain(1:end-1),chains{k}];
                        doneLookingForGreaterChain = true;
                    end
                    k = k + 1;
                end
                if ~doneLookingForGreaterChain
                    chains{s} = thisChain;
                    s = s+1;
                end
                allChains = ~any(any(S));                                
            end
        end
    end
    
    if ~isempty(chains)
        if plotPhases
            h = plot(net);
            hold on;
        end
        if length(chains) == 1 && length(chains{1}) > net.trainParam.N
            return;
        end
        % Ordering chains
        for c = 1:length(chains)
            if chains{c}(1) > chains{c}(end)
                chains{c} = fliplr(chains{c});
            end
        end

        % Consider case with unique chain and all cities included    

        % Visualize Chains
        if plotPhases
            plot(net,'phase1',chains,[],myPlot.myCitiesColorP1,myPlot.myCitiesTextColor,myPlot.myInsideColorP1,myPlot.myTitleP1);
        end
    end
end

function [netPhase2,V] = simDivideConquerPhase2(net,chains,plotPhases,myPlot)

    noChains = length(chains);   
    method = 'phase2';
    if strcmp(method, 'phase2')
        % Create new network
        % New Coordinates 
        chainExtremes  = cell(1,noChains);
        for s = 1:noChains
            chainExtremes{s}  = chains{s}([1,end]);
        end
        singleCities = setxor(1:net.trainParam.N,[chains{:}]);
        newCities = [[chainExtremes{:}],singleCities];
        
        newNames  = net.cities.names(newCities);

        newOptions = tsphopfieldnet.createOptions(...
            'R_ITER'              , net.setting.R_ITER               ,...
            'dt'                  , net.setting.dt                   ,...
            'e'                   , net.setting.e                    ,...
            'hwResources'         , net.setting.hwResources          ,...
            'loggingV'            , net.setting.loggingV             ,...
            'maxIter'             , net.setting.maxIter              ,...
            'q'                   , net.setting.q                    ,...
            'showCommandLine'     , net.setting.showCommandLine      ,...
            'u0'                  , net.setting.u0                   ,...
            'd'                   , net.cities.d(newCities,newCities),...
            'names'               , newNames                         ,...
            'type'                , net.cities.type                  ,...
            'simFcn'              , 'talavan-yanez');

        if ~isempty(net.cities.coords)
            newOptions.cities.coords = net.cities.coords(newCities,:);
        end
        
        netPhase2 = tsphopfieldnet(length(newCities),net.trainParam.C,newOptions);
        netPhase2.setting.transferFcn = net.setting.transferFcn;
        netPhase2.setting.invTransferFcn = net.setting.invTransferFcn;
        
        netPhase2.trainParam.K = noChains;
%         aux_dPhase2 = netPhase2.cities.d;
%         reOrder = reshape(flipud(reshape(1:2*netPhase2.trainParam.K,2,netPhase2.trainParam.K)),2*netPhase2.trainParam.K,1);
%         dxyc = netPhase2.cities.d;
%         dxyc(:,1:2*netPhase2.trainParam.K) = dxyc(:,reOrder);
%         dxyc(1:length(dxyc)+1:end) = 0; %?  
%         netPhase2.cities.d = dxyc;
%         netPhase2.cities.d(1,2) = 0;
%         netPhase2.cities.d(2,1) = 0;
        
        if netPhase2.trainParam.K > 1 && netPhase2.trainParam.N > 2
            for x = 1:2:2*netPhase2.trainParam.K-1
                netPhase2.cities.d(x+1,x) = 0;
                netPhase2.cities.d(x,x+1) = 0;
            end
        end
        train(netPhase2);
%         netPhase2.cities.d = aux_dPhase2;
        
        V = saddle(netPhase2) + (rand(netPhase2.trainParam.N,netPhase2.trainParam.N-netPhase2.trainParam.K) - 0.5) * 1e-10;
        U = netPhase2.setting.invTransferFcn(V);
        
        V = sim(netPhase2,V,U);
        if ~netPhase2.results.validPath
            net.results.validPath = netPhase2.results.validPath;
            net.results.exitFlag =  netPhase2.results.exitFlag;
            return;
        end
        
        if plotPhases
            fixedChainsFromPhase1 = cell(size(chains));
            for i = 1:length(fixedChainsFromPhase1)
            	fixedChainsFromPhase1{i} = [2*i-1,2*i];
            end
            plot(netPhase2,'phase2',fixedChainsFromPhase1,[],myPlot.myCitiesColorP2,myPlot.myCitiesTextColor,myPlot.myInsideColorP2,myPlot.myTitleP2);
        end

    end
    
    net.simFcn = 'divide-conquer';
    net.cities.fixedCities = {''};
    net.cities.startFixedCitiesIn = NaN;

end


% --- Auxiliar Functions for Simulation Algorithms --- %
function TV = weightMatrixTimesV(net,deltaPrev,deltaNext,...
    V, sumVrow, sumVcol, sumV, x, i)

    TV = - net.trainParam.D * net.cities.d(x,:) * ...
        sum(V(:,[deltaPrev, deltaNext]),2) + ...
        (-net.trainParam.A * sumVrow(x) - ...
        net.trainParam.B * sumVcol(i)) -...
        net.trainParam.C * sumV + ...
        (net.trainParam.A + net.trainParam.B) * V(x,i);
end

function TV = weightMatrixTimesVvectorized(net,...
    V, sumVrow, sumVcol, sumV)

    deltaPrev = [net.trainParam.N-net.trainParam.K,1:net.trainParam.N-net.trainParam.K-1];
    deltaNext = [2:net.trainParam.N-net.trainParam.K,1];
    
    if net.trainParam.K == 0
        dVpVn = net.cities.d * (V(:,deltaPrev) + V(:,deltaNext));
        TV = bsxfun(@plus, -net.trainParam.A*sumVrow, ...
            -net.trainParam.B*sumVcol) + ...
            (net.trainParam.A + net.trainParam.B) * V - ... 
            net.trainParam.C * sumV - ...
            net.trainParam.D * dVpVn;
            if isfield(net.trainParam,'F')
                TV = TV + net.trainParam.F * (net.cities.d * bsxfun(@minus, sumVrow, V) - dVpVn);
                % net.trainParam.F * net.cities.d * (repmat(sumVrow,1,net.trainParam.N) - V - V(:,deltaPrev) - V(:,deltaNext));
            end
    else % HAY QUE TRABAJAR ESTO
        reOrder = reshape(flipud(reshape(1:2*net.trainParam.K,2,net.trainParam.K)),2*net.trainParam.K,1);
    	termA = repmat(sumVrow,1,size(V,2)) - V;
        termA(1:2*net.trainParam.K,:) = termA(1:2*net.trainParam.K,:) + termA(reOrder,:);
        % Sum for free cities (x 2)
        termA(2*net.trainParam.K+1:end,:) = termA(2*net.trainParam.K+1:end,:) + termA(2*net.trainParam.K+1:end,:);
        
        termB = repmat(sumVcol,size(V,1),1) - V;

        dxyc = net.cities.d;
        dxyc(:,1:2*net.trainParam.K) = dxyc(:,reOrder);
%         dxyc(1:length(dxyc)+1:end) = 0; %?
        
        dxcy = net.cities.d;
        dxcy(1:2*net.trainParam.K,:) = dxcy(reOrder,:);
%         dxcy(1:length(dxyc)+1:end) = 0; %?
        
        termD = dxyc*V(:,deltaPrev) + dxcy*V(:,deltaNext);
        
        TV = -net.trainParam.A * termA - net.trainParam.B * termB  - ... 
            net.trainParam.C * sumV - net.trainParam.D * termD;
    end
end

function net = computeTour(net,V,iter)

    % Reconstructing original distance
    if strcmp(net.simFcn,'talavan-yanez')
        d = net.trainParam.dUaux * net.cities.d;
    else %strcmp(net.simFcn,'divide-conquer')
        d = net.cities.d;
    end
        
    V(V > 1 - power(10, -1 * net.setting.e)) = 1; %FIXME to be removed?
    V(V < power(10, -1 * net.setting.e)) = 0;
    
    if strcmp(net.simFcn,'euler') % More relaxed criteria
        V(V > 0.99) = 1;
        V(V < 0.01) = 0;
    end
    
    if all(any(V == 1, 1)) && all(any(V == 1,2)) && sum(sum(V)) == net.trainParam.N
        net.results.validPath = true;
        [~,net.results.visitOrder] = max(V);
        pos = sub2ind(size(d), net.results.visitOrder, ...
            circshift(net.results.visitOrder,[1,-1]));
        net.results.tourLength = sum(d(pos));
        net.results.exitFlag = 1;
        
    elseif net.trainParam.K > 0 && all(any(V == 1, 1)) && ...
            sum(sum(V)) == net.trainParam.N - net.trainParam.K && ...
            all([arrayfun(@(K) sum(sum(V(2*K-1:2*K,:))), 1:net.trainParam.K),sum(V(2*net.trainParam.K+1:end,:),2)'])
        
        % The computed tour has fixed chains
        net.results.validPath = true;
        [~,provisionalOrder] = max(V);
        net.results.visitOrder = zeros(1,net.trainParam.N);
        iNew = 1;
        iOld = 1;
        while iNew <= net.trainParam.N
            net.results.visitOrder(iNew) = provisionalOrder(iOld);
            if net.results.visitOrder(iNew) <= 2*net.trainParam.K
                if rem(net.results.visitOrder(iNew),2)
                    net.results.visitOrder(iNew+1) = net.results.visitOrder(iNew) + 1;
                    iNew = iNew+1;
                else
                    net.results.visitOrder(iNew+1) = net.results.visitOrder(iNew) - 1;
                    iNew = iNew+1;
                end
            end
            iNew = iNew+1;
            iOld = iOld+1;
        end
        
        pos = sub2ind(size(d), net.results.visitOrder, ...
            circshift(net.results.visitOrder,[1,-1]));
        net.results.tourLength = sum(d(pos));
        net.results.exitFlag = 1;
        
    else
        % Error produced during simulation
        if iter > net.setting.maxIter % Max iterations reached
            net.results.exitFlag = 0;
        else 
            if isempty(net.results.exitFlag)
                net.results.exitFlag = -1;
            end
        end
        net.results.validPath = false;
        net.results.tourLength = nan;                
    end
    net.results.itersReached = iter;

end
