function V = sim(net,V,U,isSaddle)
    
    if nargin == 1 % Start in the center of the Hypercube. Replace for saddle?
        U = rand(net.TrainParam.N,1)-.5; 
        V = 0.5 + 1e-7*U; 
    elseif nargin == 2
        U = net.Setting.InvTransferFcn(V);
    end
    
    % Send data to GPU
    if strcmp(net.Setting.ExecutionEnvironment,'gpu')
        U = gpuArray(U);
        V = gpuArray(V);
    end
    
    if nargin < 4
        isSaddle = false;
    end

    if ~isSaddle % Logging time and Checkpoint
        timeID = tic;

        if ~isempty(net.Setting.CheckpointPath) 
            net.Results.CheckpointFilename = ...
                utils.checkpoint.createFilename(net.Setting.CheckpointPath,net.SimFcn);
        else
            net.Results.CheckpointFilename = '';
        end
    end
    
    net = reinit(net);

    if ~isfield(net.TrainParam,'T')
        error('HopfieldNetworkGQKP:NotTrained', 'Training has not taken place yet. Use train(net).');
    end
    
    % Classic Scheme support only
    if strcmp(net.SimFcn,'euler')
        [net,V,~,iter] = simEuler(net,V,U,isSaddle);

    elseif strcmp(net.SimFcn,'talavan-yanez') 
        [net,V,~,iter] = simTalavanYanez(net,V,U,isSaddle);
    end
    net = computeSolution(net,V,iter);

    net.Results.CompTime = toc(timeID);
end


% Euler Algorithm
function [net,V,U,iter] = simEuler(net,V,U,isSaddle)

    % Stopping criteria
    stopC1 = power(10, -1   * net.Setting.E);
    stopC2 = power(10, -1.5 * net.Setting.E);
    maxDiffV = 1;
    unstable = false;

    N = net.ProblemParameters.networkSize(1);
    T = net.TrainParam.T;
    ib = net.TrainParam.ib;

    % Memory allocation
    dU = zeros(N,1);

    iter = 1;
% Showing Network parameters in the command window
%             if net.Setting.Verbose
%                 hopfield.tsp.display.printWhenStarting(net);
%             end

    if ~isSaddle % Logging Checkpoint and plotting Simulation process
        if ~isempty(net.Setting.CheckpointPath) || net.Setting.SimulationPlot
            if ~isempty(net.Setting.CheckpointPath) 
                utils.checkpoint.loggingData(fullfile(net.Setting.CheckpointPath, ...
                    net.Results.CheckpointFilename),...
                    net.Setting.MaxIter,iter,V,zeros(size(V)));
            end
            if net.Setting.SimulationPlot
                fV = viewConvergence(iter,V,net);
            end
        end
    end
    
    dt = net.Setting.Dt;
    net.Results.Time(iter) = dt;
    maxDiffV = 1;

    while iter <= net.Setting.MaxIter && (maxDiffV > stopC1 || ...
            (maxDiffV > stopC2 && unstable))

        unstable = false;
        net.Results.Energy(iter+1) = 0;

        dU = T*V + ib;

        Vprev = V;
        U = U + dU * net.Setting.Dt;
        V = net.Setting.TransferFcn(U);

        maxDiffV = max(abs(Vprev-V));
        iter = iter + 1;
        
        if ~isSaddle % Logging Checkpoint and plotting Simulation process
            if ~isempty(net.Setting.CheckpointPath) || net.Setting.SimulationPlot
                if ~isempty(net.Setting.CheckpointPath) 
                    utils.checkpoint.loggingData(fullfile(net.Setting.CheckpointPath,...
                        net.Results.CheckpointFilename),...
                        net.Setting.MaxIter,iter,V,dU);
                end
                if net.Setting.SimulationPlot
                    fV = viewConvergence(iter,V,net,fV);
                end
            end
        end        
    end
    
    if strcmp(net.Setting.ExecutionEnvironment,'gpu') && nargout > 1
        V = gather(V);
        U = gather(U);
    end

    if ~isSaddle % Logging Checkpoint 
        if ~isempty(net.Setting.CheckpointPath)
            utils.checkpoint.trimToSimulatedData(net.Setting.CheckpointPath, ...
                net.Results.CheckpointFilename, iter)
        end
    end
    
    net.Results.x = V;
end

% Talaván-Yáñez Algorithm

function [net,V,U,iter] = simTalavanYanez(net,V,U,isSaddle)

    T = net.TrainParam.T;
    ib = net.TrainParam.ib;
       
    % Stopping criteria
    stopC1 = power(10, -1   * net.Setting.E);
    stopC2 = power(10, -1.5 * net.Setting.E);
    maxDiffV = 1;
    unstable = false;
    
    trasferFcn2Str = func2str(net.Setting.TransferFcn);
    if strcmp(trasferFcn2Str,'@(u)0.5*(1+tanh(u./net.Setting.U0))')
        trasferFcn2Str = 'tanh';
    elseif strcmp(trasferFcn2Str,'@(u)net.satlin(u,net.Setting.U0)')
        trasferFcn2Str = 'satlin';
    else
        error('HopfieldNetworkTSP:sim:notDefinedTransferFcn', ...
            ['TransferFcn derivative not defined for ',net.Setting.TransferFcn]);
    end    
    
    u_e = net.Setting.InvTransferFcn(stopC1); 

    iter = 1;
% Showing Network parameters in the command window
%             if net.Setting.Verbose
%                 hopfield.tsp.display.printWhenStarting(net);
%             end

    if ~isSaddle % Logging Checkpoint and plotting Simulation process
        if ~isempty(net.Setting.CheckpointPath) || net.Setting.SimulationPlot
            if ~isempty(net.Setting.CheckpointPath) 
                utils.checkpoint.loggingData(fullfile(net.Setting.CheckpointPath, ...
                    net.Results.CheckpointFilename),...
                    net.Setting.MaxIter,iter,V,zeros(size(V)));
            end
            if net.Setting.SimulationPlot
                fV = viewConvergence(iter,V,net);
            end
        end
    end

    while iter < net.Setting.MaxIter && (maxDiffV > stopC1 || ...
            (maxDiffV > stopC2 && unstable))

        dt = 10^100; % Initial value for dt

        % Computation of the weight matrix T
        % $$ T_{xi,yj} = -(A*\delta_{x,y}*(1-\delta_{i,j}) + B*(1-\delta_{x,y})*
        % \delta_{i,j} + C - D*d_{x*y} * (\delta_{j,i+1} + \delta_{j,i-1}) $$

        dU = T*V + ib;

        if strcmp(trasferFcn2Str,'tanh')
            dV = 2./net.Setting.U0 .* V .* (1-V) .*dU;
        elseif strcmp(trasferFcn2Str,'satlin')
            dV = 2./net.Setting.U0 .* dU;
        end
        
        interiorV = U > u_e & U < -u_e;

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

        TdV = T*dV;

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

        %%%
        % $E(t + \Delta t) = E(t) - S_{1}\Delta t + \frac{1}{2}S_{2}\Delta t^{2}$
        if strcmp(net.Setting.ExecutionEnvironment,'gpu')
            S1 = gather(S1);
            S2 = gather(S2);
            dt = gather(dt);
        end
        net.Results.Energy(iter+1) = net.Results.Energy(iter) - ...
            S1.*dt + 0.5.*S2.*dt.^2;
        net.Results.Time(iter+1) = net.Results.Time(iter) + dt;
        iter = iter + 1;
        
        if ~isSaddle % Logging Checkpoint and plotting Simulation process
            if ~isempty(net.Setting.CheckpointPath) || net.Setting.SimulationPlot
                if ~isempty(net.Setting.CheckpointPath) 
                    utils.checkpoint.loggingData(fullfile(net.Setting.CheckpointPath,...
                        net.Results.CheckpointFilename),...
                        net.Setting.MaxIter,iter,V,dU);
                end
                if net.Setting.SimulationPlot
                    fV = viewConvergence(iter,V,net,fV);
                end
            end
        end
    end

    if strcmp(net.Setting.ExecutionEnvironment,'gpu') && nargout > 1
        V = gather(V);
        U = gather(U);
    end

    if ~isSaddle % Logging Checkpoint 
        if ~isempty(net.Setting.CheckpointPath)
            utils.checkpoint.trimToSimulatedData(net.Setting.CheckpointPath, ...
                net.Results.CheckpointFilename, iter);
        end
    end   
    
end

function net = computeSolution(net,V,iter)

    V(V > 1 - power(10, -1 * net.Setting.E)) = 1; 
    V(V < power(10, -1 * net.Setting.E)) = 0;

    if isequal(net.ProblemParameters.R*V,net.ProblemParameters.b)
        net.Results.ValidSolution = true;
        net.Results.ExitFlag = 1;
    else
        % Error produced during simulation
        if iter > net.Setting.MaxIter % Max iterations reached
            net.Results.ExitFlag = 0;
        else
            net.Results.ExitFlag = -1;
        end
        net.Results.ValidSolution = false;
    end
    net.Results.x = V;
	net.Results.fval = 0.5*net.Results.x'*net.ProblemParameters.P*net.Results.x + net.ProblemParameters.q'*net.Results.x;

    net.Results.ItersReached = iter;

end
