function V = sim(net,V,U)

	[ST,~] = dbstack('-completenames');
    if ~any(strcmp({ST.name},'saddle'))
        timeID = tic;

        if ~isempty(net.Setting.CheckpointPath) 
            warning('off','MATLAB:MKDIR:DirectoryExists');
            mkdir(net.Setting.CheckpointPath);
            warning('on','MATLAB:MKDIR:DirectoryExists');
            checkpointFilename = ['checkpoint_',net.SimFcn,'_',...
                datestr(datetime(now,'ConvertFrom','datenum'),...
                        'yyyy_mm_dd_HH_MM_SS'),'.mat'];
            net.Results.CheckpointFilename = checkpointFilename;
        else
            net.Results.CheckpointFilename = '';
        end
    end
    
    net = reinit(net);
    
    if ~isfield(net.TrainParam,'Np')
        error('tsphopfieldnet:NotTrained', 'Training has not taken place yet. Use train(net).');
    end
       
    if ~isempty(net.Cities.Subtours)
        tsphopfieldnet.verifyIfValidSubtours(net.Cities.Subtours, net.Cities.SubtoursPositions, net.Cities.Names)
        aux_d = net.Cities.DistanceMatrix;
        [net,V,U] = fixedCities(net); %Review fixedCities
    end
    
    if strcmp(net.SimFcn,'euler')
        if nargin < 2
            [net,V,~,iter] = simEuler(net);
        else
            [net,V,~,iter] = simEuler(net,V,U);
        end
        
    elseif strcmp(net.SimFcn,'talavan-yanez')
        if nargin < 2 && isempty(net.Cities.Subtours)
            [net,V,~,iter] = simTalavanYanez(net);
        else
            [net,V,~,iter] = simTalavanYanez(net,V,U);
        end
        
    elseif strcmp(net.SimFcn, 'divide-conquer')
        if nargin < 2
            [net,V,~,iter] = simDivideConquer(net);
        else
            [net,V,~,iter] = simDivideConquer(net,V,U);
        end
        
    else
        error('tsphopfieldnet:UnknownsimFcn', 'Unknown training algorithm');
    end
    
    if ~isempty(net.Cities.Subtours)
        net.Cities.DistanceMatrix = aux_d;
    end   
    
    if ~any(strcmp({ST.name},'saddle'))
        net = computeTour(net,V,iter);
        net.Results.CompTime = toc(timeID);

        % Removing unused energy and time elements.
        net.Results.Energy = net.Results.Energy(~isnan(net.Results.Energy));
        indexRemove = find(diff(net.Results.Time),1,'last') + 1;
        if length(net.Results.Time) > indexRemove
            net.Results.Time(indexRemove:end) = [];
        end
%         net.Results.time = [net.Results.time(1),net.Results.time(net.Results.time(2:end) ~= 0)];

        if ~isempty(net.Setting.CheckpointPath)
            data = matfile(fullfile(net.Setting.CheckpointPath,...
            net.Results.CheckpointFilename),'Writable',true);
            data.net = net;
        end
        net.Results = orderfields(net.Results);
    end
end

% --- Simulation Algorithms --- %

% Euler Algorithm
function [net,V,U,iter] = simEuler(net, V, U) 
    N = net.TrainParam.N;
    if nargin == 1
        U = rand(N)-.5;   % TODO Different from Pedro's algorithm
        V = 0.5 + 1e-7*U; % TODO Different from Pedro's algorithm
    end
    if strcmp(net.Setting.ExecutionEnvironment,'GPU')
        U = gpuArray(U);
        V = gpuArray(V);
    end
    % Stopping criteria
    stopC1 = power(10, -1   * net.Setting.E);
    stopC2 = power(10, -1.5 * net.Setting.E);
    maxDiffV = 1;
    unstable = false;

    ib = net.TrainParam.C * net.TrainParam.Np;

    %%%
    % Summing V rows, columns and all elements for optimal perform.
    sumVcol = sum(V);
    sumVrow = sum(V,2);
    sumV = sum(sumVcol);

    % Memory allocation
    dU = zeros(N);

    iter = 1;

    dt = net.Setting.Dt;
    net.Results.Time(iter) = dt;

    if ~isempty(net.Setting.CheckpointPath) || net.Setting.SimulationPlot
        if ~isempty(net.Setting.CheckpointPath) 
            tsphopfieldnet.loggingV(fullfile(net.Setting.CheckpointPath, ...
                net.Results.CheckpointFilename),...
                net.Setting.MaxIter,iter,N);
        end
        if net.Setting.SimulationPlot
            fV = viewConvergence(iter,V,net);
        end
    end    
    
    while iter <= net.Setting.MaxIter && (maxDiffV > stopC1 || ...
            (maxDiffV > stopC2 && unstable))

        unstable = false;
        for i = 1:N
            if i == N, deltaPrev = 1; else, deltaPrev = i+1; end
            if i == 1, deltaNext = N; else, deltaNext = i-1; end
            for x = 1:N
                TV = weightMatrixTimesV(...
                    net, deltaPrev, deltaNext, V, sumVrow, ...
                    sumVcol, sumV, x, i);

                dU(x,i) = TV + ib;

                if (V(x,i) < stopC1 && dU(x,i) > 0) || ...
                        (V(x,i) > 1 - stopC1 && dU(x,i) < 0)
                    unstable = true;
                end
                net.Results.Energy(iter+1) = ...
                    net.Results.Energy(iter) + ...
                    0.5* V(x,i) * dU(x,i) - ...
                    ib * V(x,i);
            end
        end

        maxDiffV = 0;
        for i = 1:N
            for x = 1:N
                prevV = V(x,i);
                U(x,i) = U(x,i) + dU(x,i)*net.Setting.Dt;
                V(x,i) = net.Setting.TransferFcn(U(x,i));
                if abs(prevV - V(x,i)) > maxDiffV
                    maxDiffV = abs(prevV - V(x,i));
                end
            end          
        end   
        sumVcol = sum(V);
        sumVrow = sum(V,2);
        sumV = sum(sumVcol);       

        net.Results.Time(iter+1) = net.Results.Time(iter) + ...
            net.Setting.Dt;
        iter = iter + 1;
        if ~isempty(net.Setting.CheckpointPath) || net.Setting.SimulationPlot
            if ~isempty(net.Setting.CheckpointPath) 
                tsphopfieldnet.loggingV(fullfile(net.Setting.CheckpointPath,...
                    net.Results.CheckpointFilename),...
                    net.Setting.MaxIter,iter,[],V,dU);
            end
            if net.Setting.SimulationPlot
                fV = viewConvergence(iter,V,net,fV);
            end
        end        
    end
    if ~isempty(net.Setting.CheckpointPath)
        data = matfile(fullfile(net.Setting.CheckpointPath,...
            net.Results.CheckpointFilename),'Writable',true);
        data.Vlog(:,:,iter+1:end) = [];
        data.dUlog(:,:,iter+1:end) = [];
    end    

end

% Talaván-Yáñez Algorithm
function [net,V,U,iter] = simTalavanYanez(net,V,U)
    N = net.TrainParam.N;
    K = net.TrainParam.K;
    if nargin == 1
        U = rand(N,N-K)-.5;   % TODO Different from Pedro's algorithm
        V = 0.5 + 1e-7*U;     % TODO Different from Pedro's algorithm
    end
    if strcmp(net.Setting.ExecutionEnvironment,'GPU')
        U = gpuArray(U);
        V = gpuArray(V);
    end
    % Stopping criteria
    stopC1 = power(10, -1   * net.Setting.E);
    stopC2 = power(10, -1.5 * net.Setting.E);
    maxDiffV = 1;
    unstable = false;

    u_e = net.Setting.InvTransferFcn(stopC1); 

	ib = net.TrainParam.C * net.TrainParam.Np;
    
    %%%
    % Summing V rows, columns and all elements for optimal performance
    sumVcol = sum(V);   % Si
    sumVrow = sum(V,2); % Sx
    sumV = sum(sumVcol);

    % Memory allocation
    dU = zeros(N,N-K);

    iter = 1;
 
    [ST,~] = dbstack('-completenames');
    if ~any(strcmp({ST.name},'saddle'))
        if ~isempty(net.Setting.CheckpointPath) || net.Setting.SimulationPlot
            if ~isempty(net.Setting.CheckpointPath) 
                tsphopfieldnet.loggingV(fullfile(net.Setting.CheckpointPath, ...
                    net.Results.CheckpointFilename),...
                    net.Setting.MaxIter,iter,N-K);
            end
            if net.Setting.SimulationPlot
                fV = viewConvergence(iter,V,net);
            end
        end
    end
           
    while iter <= net.Setting.MaxIter && (maxDiffV > stopC1 || ...
            (maxDiffV > stopC2 && unstable))

        dt = 10^100; % Initial value for dt

        % Computation of the weight matrix T
        % $$ T_{xi,yj} = -(A*\delta_{x,y}*(1-\delta_{i,j}) + B*(1-\delta_{x,y})*
        % \delta_{i,j} + C - D*d_{x*y} * (\delta_{j,i+1} + \delta_{j,i-1}) $$
        TV = weightMatrixTimesVvectorized(net, V, sumVrow, sumVcol, sumV);
        
        dU = TV + ib;

        dV = 2./net.Setting.U0 .* V .* (1-V) .*dU;

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
        if iter < net.Setting.R_Iter && ~sw_optimal
            dt = net.Setting.Q * dt; % Integration step reduction
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
        V(borderV) = net.Setting.TransferFcn(U(borderV));
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
        if strcmp(net.Setting.ExecutionEnvironment,'GPU')
            S1 = gather(S1);
            S2 = gather(S2);
            dt = gather(dt);
        end
        net.Results.Energy(iter+1) = net.Results.Energy(iter) - ...
            S1.*dt + 0.5.*S2.*dt.^2;
        net.Results.Time(iter+1) = net.Results.Time(iter) + dt;
        iter = iter + 1;
        
        if ~any(strcmp({ST.name},'saddle'))
            if ~isempty(net.Setting.CheckpointPath) || net.Setting.SimulationPlot
                if ~isempty(net.Setting.CheckpointPath) 
                    tsphopfieldnet.loggingV(fullfile(net.Setting.CheckpointPath,...
                        net.Results.CheckpointFilename),...
                        net.Setting.maxIter,iter,[],V,dU);
                end
                if net.Setting.SimulationPlot
                    fV = viewConvergence(iter,V,net,fV);
                end
            end
        end
    end

    if strcmp(net.Setting.ExecutionEnvironment,'GPU') && nargout > 1
        V = gather(V);
        U = gather(U);
        iter = gather(iter);
    end

    if ~any(strcmp({ST.name},'saddle'))
        if ~isempty(net.Setting.CheckpointPath)
            data = matfile(fullfile(net.Setting.CheckpointPath,...
                net.Results.CheckpointFilename),'Writable',true);
            data.Vlog(:,:,iter+1:end) = [];
            data.dUlog(:,:,iter+1:end) = [];
        end
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
    %    b. using p which uses distances <= p * net.TrainParam.dU + (1-p) * net.TrainParam.dL
    
    if net.TrainParam.K == 0
        if nargin < 2
            V = saddle(net) + (rand(net.TrainParam.N) - 0.5) * 1e-5;
            U = net.Setting.InvTransferFcn(V);
        end
        [chains,V] = simDivideConquerPhase1(net,V,U,plotPhases,myPlot);
        
        % Reached max number of iterations        
        if net.Results.ExitFlag == 0
            mkdir('Phase1_InsuficientIterations');
            save(fullfile(pwd,'Phase1_InsuficientIterations',regexprep(datestr(datetime),{':',' ','-'},'_')))
            while net.Results.ExitFlag <= 0
                V = V + (rand(net.TrainParam.N) - 0.5) * 1e-5;
                U = net.Setting.InvTransferFcn(V);
                [chains,V] = simDivideConquerPhase1(net,V,U,plotPhases,myPlot);
            end
            
        % Possible saddle point reached
        elseif net.Results.ExitFlag == -1 
            mkdir('Phase1_PossibleSaddlePoint');
            save(fullfile(pwd,'Phase1_PossibleSaddlePoint',regexprep(datestr(datetime),{':',' ','-'},'_')))            
            while net.Results.ExitFlag <= 0
                V = V + (rand(net.TrainParam.N) - 0.5) * 1e-5;
                U = net.Setting.InvTransferFcn(V);
                [chains,V] = simDivideConquerPhase1(net,V,U,plotPhases,myPlot);
            end
        end
        % TODO Add correct number of iterations
    else
        % Construir chains
        chains = cell(1,net.TrainParam.K);
        for c = 1:net.TrainParam.K
            chains{c} = [2*c-1,2*c];
        end
    end
    
    if (~isempty(chains) && length(chains{1}) <= net.TrainParam.N && net.Results.ExitFlag == 1) || isempty(chains) && net.Results.ExitFlag == 1

        % Part 2. 3 possible methods:
        % a. fixing cities
        % b. using new distance
        % c. connecting subtours with greedy
        [netPhase2,V2] = simDivideConquerPhase2(net,chains,plotPhases,myPlot);
        % Building final V
        V = zeros(net.TrainParam.N);
        iNew = 1;
        iOld = 1;

        % Reached max number of iterations        
        if netPhase2.Results.ExitFlag == 0
            mkdir('Phase2_InsuficientIterations');
            save(fullfile(pwd,'Phase2_InsuficientIterations',regexprep(datestr(datetime),{':',' ','-'},'_')))
            while netPhase2.Results.ExitFlag <= 0
                V2 = V2 + (rand(netPhase2.TrainParam.N,netPhase2.TrainParam.N-netPhase2.TrainParam.K) - 0.5) * 1e-5;
                U = netPhase2.Setting.InvTransferFcn(V2);
                V2 = sim(netPhase2, V2, U);
            end
            
        % Possible saddle point reached
        elseif netPhase2.Results.ExitFlag == -1 
            mkdir('Phase2_PossibleSaddlePoint');
            save(fullfile(pwd,'Phase2_PossibleSaddlePoint',regexprep(datestr(datetime),{':',' ','-'},'_')))            
            % Stuck in a saddle point. Perturbation required.
            while netPhase2.Results.ExitFlag <= 0
                V2 = V2 + (rand(netPhase2.TrainParam.N,netPhase2.TrainParam.N-netPhase2.TrainParam.K) - 0.5) * 1e-5;
                U = netPhase2.Setting.InvTransferFcn(V2);
                V2 = sim(netPhase2, V2, U);
            end
        end
        % TODO Add correct number of iterations        
        
        while iNew <= net.TrainParam.N
            
            thisCity = find(strcmp(net.Cities.Names,netPhase2.Cities.Names(netPhase2.Results.VisitOrder(iOld))));

            V(thisCity,iNew) = 1;
            iNew = iNew + 1;
            if netPhase2.Results.VisitOrder(iOld) <= 2*netPhase2.TrainParam.K
                for c = 1:length(chains)
                    if any(chains{c} == thisCity)
                        break;
                    end
                end
                if rem(netPhase2.Results.VisitOrder(iOld),2) % Start of chain
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
        U = net.Setting.InvTransferFcn(V);
        iter = net.Results.ItersReached + netPhase2.Results.ItersReached - 1;

        if plotPhases
            if isempty(net.Results.TourLength) %#ok<UNRCH>
                init(net);
                net = computeTour(net,V,iter); % Needed for plot phase1 + phase2 to output correctly
            else
                [~,net.Results.VisitOrder] = max(V); % Needed for plot phase1 + phase2 to output correctly
            end
            plot(net,'phase2',chains,[],myPlot.myCitiesColorP2,myPlot.myCitiesTextColor,myPlot.myInsideColorP1,myPlot.myTitle);
        end    

        % Removing unused energy and time elements.
        net.Results.Energy = net.Results.Energy(~isnan(net.Results.Energy));
        net.Results.Time = [net.Results.Time(1),net.Results.Time(net.Results.Time(2:end) ~= 0)];

        net.Results.Energy = [net.Results.Energy, netPhase2.Results.Energy(2:end)];
        net.Results.Time = [net.Results.Time, net.Results.Time(end) + netPhase2.Results.Time(2:end)];
    else
        iter = net.Results.ItersReached;
    end
end

function [chains,V] = simDivideConquerPhase1(net,V,U,plotPhases,myPlot)

    p_or_tau = net.Cities.Tau;
	% Backing up distances
    aux_d = net.Cities.DistanceMatrix;

    [net.Cities.d, neighbours] = neighbourDistance(net, p_or_tau);

    % Starting point close to saddle point
    if nargin < 2
        V = saddle(net) + (rand(net.TrainParam.N) - 0.5) * 1e-5;
        U = net.Setting.InvTransferFcn(V);
    end
    
    [net,V,U,iter] = simTalavanYanez(net,V,U);

    % Retrieving original distance from backup.    
    net.Cities.d = aux_d;
    
    % Finding chains such that go from city x to city y (from
    % position j to position j+1) using the original distance.
   
    net = computeTour(net,V,iter);
    % Extract chains from S.
	chains = {};
            
    if net.Results.ValidPath
        finalTour = net.Results.VisitOrder;
        
        % X: City matrix. Connecting cities from j to j+1
        X = zeros(net.TrainParam.N);
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
        if length(chains) == 1 && length(chains{1}) > net.TrainParam.N
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
        singleCities = setxor(1:net.TrainParam.N,[chains{:}]);
        newCities = [[chainExtremes{:}],singleCities];
        
        newNames  = net.Cities.Names(newCities);
        if ~isempty(net.Cities.Coordinates)
            newCoordinates = net.Cities.Coordinates(newCities,:);
        else
            newCoordinates = [];
        end
        
        newOptions = tsphopfieldnetOptions(...
            'R_Iter'              , net.Setting.R_Iter                            ,...
            'Dt'                  , net.Setting.Dt                                ,...
            'E'                   , net.Setting.E                                 ,...
            'ExecutionEnvironment', net.Setting.ExecutionEnvironment              ,...
            'CheckpointPath'      , net.Setting.CheckpointPath                    ,...
            'MaxIter'             , net.Setting.MaxIter                           ,...
            'Q'                   , net.Setting.Q                                 ,...
            'Verbose'             , net.Setting.Verbose                           ,...
            'U0'                  , net.Setting.U0                                ,...
            'DistanceMatrix'      , net.Cities.DistanceMatrix(newCities,newCities),...
            'Names'               , newNames                                      ,...
            'DistanceType'        , net.Cities.DistanceType                       ,...
            'Coordinates'         , newCoordinates                                ,...
            'SimFcn'              , 'talavan-yanez');


        
        netPhase2 = tsphopfieldnet(length(newCities),net.TrainParam.C,newOptions);
        netPhase2.Setting.transferFcn = net.Setting.TransferFcn;
        netPhase2.Setting.InvTransferFcn = net.Setting.InvTransferFcn;
        
        netPhase2.TrainParam.K = noChains;
%         aux_dPhase2 = netPhase2.cities.d;
%         reOrder = reshape(flipud(reshape(1:2*netPhase2.TrainParam.K,2,netPhase2.TrainParam.K)),2*netPhase2.TrainParam.K,1);
%         dxyc = netPhase2.cities.d;
%         dxyc(:,1:2*netPhase2.TrainParam.K) = dxyc(:,reOrder);
%         dxyc(1:length(dxyc)+1:end) = 0; %?  
%         netPhase2.cities.d = dxyc;
%         netPhase2.cities.d(1,2) = 0;
%         netPhase2.cities.d(2,1) = 0;
        
        if netPhase2.TrainParam.K > 1 && netPhase2.TrainParam.N > 2
            for x = 1:2:2*netPhase2.TrainParam.K-1
                netPhase2.Cities.DistanceMatrix(x+1,x) = 0;
                netPhase2.Cities.DistanceMatrix(x,x+1) = 0;
            end
        end
        train(netPhase2);
%         netPhase2.cities.d = aux_dPhase2;
        
        V = saddle(netPhase2) + (rand(netPhase2.TrainParam.N,netPhase2.TrainParam.N-netPhase2.TrainParam.K) - 0.5) * 1e-10;
        U = netPhase2.Setting.InvTransferFcn(V);
        
        V = sim(netPhase2,V,U);
        if ~netPhase2.Results.ValidPath
            net.Results.ValidPath = netPhase2.Results.ValidPath;
            net.Results.ExitFlag =  netPhase2.Results.ExitFlag;
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
    
    net.SimFcn = 'divide-conquer';
    net.Cities.Subtours = '';
    net.Cities.SubtoursPositions = NaN;

end


% --- Auxiliar Functions for Simulation Algorithms --- %
function TV = weightMatrixTimesV(net,deltaPrev,deltaNext,...
    V, sumVrow, sumVcol, sumV, x, i)

    TV = - net.TrainParam.D * net.Cities.DistanceMatrix(x,:) * ...
        sum(V(:,[deltaPrev, deltaNext]),2) + ...
        (-net.TrainParam.A * sumVrow(x) - ...
        net.TrainParam.B * sumVcol(i)) -...
        net.TrainParam.C * sumV + ...
        (net.TrainParam.A + net.TrainParam.B) * V(x,i);
end

function TV = weightMatrixTimesVvectorized(net,...
    V, sumVrow, sumVcol, sumV)

    deltaPrev = [net.TrainParam.N-net.TrainParam.K,1:net.TrainParam.N-net.TrainParam.K-1];
    deltaNext = [2:net.TrainParam.N-net.TrainParam.K,1];
    
    if net.TrainParam.K == 0
        dVpVn = net.Cities.DistanceMatrix * (V(:,deltaPrev) + V(:,deltaNext));
        TV = bsxfun(@plus, -net.TrainParam.A*sumVrow, ...
            -net.TrainParam.B*sumVcol) + ...
            (net.TrainParam.A + net.TrainParam.B) * V - ... 
            net.TrainParam.C * sumV - ...
            net.TrainParam.D * dVpVn;

    else % HAY QUE TRABAJAR ESTO
        reOrder = reshape(flipud(reshape(1:2*net.TrainParam.K,2,net.TrainParam.K)),2*net.TrainParam.K,1);
    	termA = repmat(sumVrow,1,size(V,2)) - V;
        termA(1:2*net.TrainParam.K,:) = termA(1:2*net.TrainParam.K,:) + termA(reOrder,:);
        % Sum for free cities (x 2)
        termA(2*net.TrainParam.K+1:end,:) = termA(2*net.TrainParam.K+1:end,:) + termA(2*net.TrainParam.K+1:end,:);
        
        termB = repmat(sumVcol,size(V,1),1) - V;

        dxyc = net.Cities.DistanceMatrix;
        dxyc(:,1:2*net.TrainParam.K) = dxyc(:,reOrder);
%         dxyc(1:length(dxyc)+1:end) = 0; %?
        
        dxcy = net.Cities.DistanceMatrix;
        dxcy(1:2*net.TrainParam.K,:) = dxcy(reOrder,:);
%         dxcy(1:length(dxyc)+1:end) = 0; %?
        
        termD = dxyc*V(:,deltaPrev) + dxcy*V(:,deltaNext);
        
        TV = -net.TrainParam.A * termA - net.TrainParam.B * termB  - ... 
            net.TrainParam.C * sumV - net.TrainParam.D * termD;
    end
end

function net = computeTour(net,V,iter)

    % Reconstructing original distance
    if strcmp(net.SimFcn,'talavan-yanez')
        d = net.TrainParam.dUaux * net.Cities.DistanceMatrix;
    else %strcmp(net.SimFcn,'divide-conquer')
        d = net.Cities.DistanceMatrix;
    end
        
    V(V > 1 - power(10, -1 * net.Setting.E)) = 1; %FIXME to be removed?
    V(V < power(10, -1 * net.Setting.E)) = 0;
    
    if strcmp(net.SimFcn,'euler') % More relaxed criteria
        V(V > 0.99) = 1;
        V(V < 0.01) = 0;
    end
    
    if all(any(V == 1, 1)) && all(any(V == 1,2)) && sum(sum(V)) == net.TrainParam.N
        net.Results.ValidPath = true;
        [~,net.Results.VisitOrder] = max(V);
        pos = sub2ind(size(d), net.Results.VisitOrder, ...
            circshift(net.Results.VisitOrder,[1,-1]));
        net.Results.TourLength = sum(d(pos));
        net.Results.ExitFlag = 1;
        
    elseif net.TrainParam.K > 0 && all(any(V == 1, 1)) && ...
            sum(sum(V)) == net.TrainParam.N - net.TrainParam.K && ...
            all([arrayfun(@(K) sum(sum(V(2*K-1:2*K,:))), 1:net.TrainParam.K),sum(V(2*net.TrainParam.K+1:end,:),2)'])
        
        % The computed tour has fixed chains
        net.Results.ValidPath = true;
        [~,provisionalOrder] = max(V);
        net.Results.VisitOrder = zeros(1,net.TrainParam.N);
        iNew = 1;
        iOld = 1;
        while iNew <= net.TrainParam.N
            net.Results.VisitOrder(iNew) = provisionalOrder(iOld);
            if net.Results.VisitOrder(iNew) <= 2*net.TrainParam.K
                if rem(net.Results.VisitOrder(iNew),2)
                    net.Results.VisitOrder(iNew+1) = net.Results.VisitOrder(iNew) + 1;
                    iNew = iNew+1;
                else
                    net.Results.VisitOrder(iNew+1) = net.Results.VisitOrder(iNew) - 1;
                    iNew = iNew+1;
                end
            end
            iNew = iNew+1;
            iOld = iOld+1;
        end
        
        pos = sub2ind(size(d), net.Results.VisitOrder, ...
            circshift(net.Results.VisitOrder,[1,-1]));
        net.Results.TourLength = sum(d(pos));
        net.Results.ExitFlag = 1;
        
    else
        % Error produced during simulation
        if iter > net.Setting.MaxIter % Max iterations reached
            net.Results.ExitFlag = 0;
        else 
            if isempty(net.Results.ExitFlag)
                net.Results.ExitFlag = -1;
            end
        end
        net.Results.ValidPath = false;
        net.Results.TourLength = nan;                
    end
    net.Results.ItersReached = iter;

end
