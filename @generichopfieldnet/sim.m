function V = sim(net,V,U)

    timeID = tic;
    net = reinit(net);
    if strcmp(net.simFcn,'euler')
        if nargin < 2
            [net,V,~,iter] = simEuler(net);
        else
            [net,V,~,iter] = simEuler(net,V,U);
        end  
    elseif strcmp(net.simFcn,'talavan-yanez') 
        if nargin < 2
            [net,V,~,iter] = simTalavanYanez(net);
        else
            [net,V,~,iter] = simTalavanYanez(net,V,U);
        end
    end
    net = computeSolution(net,V,iter);

    net.results.compTime = toc(timeID);
end


% Euler Algorithm

function [net,V,U,iter] = simEuler(net,V)

    % Stopping criteria
    stopC1 = power(10, -1   * net.setting.e);
    stopC2 = power(10, -1.5 * net.setting.e);
    maxDiffV = 1;
    unstable = false;

    %%%
    % Summing V rows, columns and all elements for optimal perform.
%     sumVcol = sum(V);
%     sumVrow = sum(V,2);
%     sumV = sum(sumVcol);

    N = net.problemParameters.networkSize(1);
    T = net.trainParam.T;
    ib = net.trainParam.ib;
    net.setting.maxIter = 10000;
    % Memory allocation
    dU = zeros(N,1);

    if nargin == 1
        U = rand(N,1)-.5;   % TODO Different from Pedro's algorithm
        V = 0.5 + 1e-7*U; % TODO Different from Pedro's algorithm
    end

    iter = 1;
    % Showing Network parameters in the command window
%             if net.setting.showCommandLine
%                 hopfield.tsp.display.printWhenStarting(net);
%             end

    dt = net.setting.dt;
    net.results.time(iter) = dt;
    maxDiffV = 1;

    while iter <= net.setting.maxIter && (maxDiffV > stopC1 || ...
            (maxDiffV > stopC2 && unstable))

        unstable = false;
        net.results.energy(iter+1) = 0;

        dU = T*V + ib;

        Vprev = V;
        U = U + dU * net.setting.dt;
        V = net.setting.transferFcn(U);

        maxDiffV = max(abs(Vprev-V));
        iter = iter + 1;            
    end
    net.results.x = V;
end

% Talaván-Yáñez Algorithm

function [net,V,U,iter] = simTalavanYanez(net,V,U)
%     N = net.trainParam.N;

    N = net.problemParameters.networkSize(1);
    T = net.trainParam.T;
    ib = net.trainParam.ib;
    
    if nargin == 1
        U = rand(N,1)-.5;   % TODO Different from Pedro's algorithm
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

    u_e = net.setting.invTransferFcn(stopC1); 

%     ib = net.trainParam.C * net.trainParam.Np;

    %%%
    % Summing V rows, columns and all elements for optimal performance
%     sumVcol = sum(V);   % Si
%     sumVrow = sum(V,2); % Sx
%     sumV = sum(sumVcol);

    % Memory allocation
%             dU = zeros(N);

    iter = 1;
    % Showing Network parameters in the command window
%             if net.setting.showCommandLine
%                 hopfield.tsp.display.printWhenStarting(net);
%             end

    while iter < net.setting.maxIter && (maxDiffV > stopC1 || ...
            (maxDiffV > stopC2 && unstable))

        dt = 10^100; % Initial value for dt

%         if net.setting.loggingV || net.setting.viewConvergence
%             [ST,~] = dbstack('-completenames');
%             if ~any(strcmp({ST.name},'saddle'))
%                 if net.setting.loggingV
%                     tsphopfieldnet.loggingV(iter,V);
%                 end
%                 if net.setting.viewConvergence
%                     if iter == 1
%                         fV = viewConvergence(iter,V,net);
%                     else
%                         fV = viewConvergence(iter,V,net,fV);
%                     end
%                 end
%             end
%         end

        % Computation of the weight matrix T
        % $$ T_{xi,yj} = -(A*\delta_{x,y}*(1-\delta_{i,j}) + B*(1-\delta_{x,y})*
        % \delta_{i,j} + C - D*d_{x*y} * (\delta_{j,i+1} + \delta_{j,i-1}) $$
%         TV = weightMatrixTimesVvectorized(net, V, ...
%             sumVrow, sumVcol, sumV);

        dU = T*V + ib;

        dV = 2./net.setting.u0 .* V .* (1-V) .*dU;

%                 interiorV = V > 0 & V < 1;
        interiorV = U > u_e & U < -u_e;
        %(V > stopC1) & (V < 1-stopC1);
%                 borderV = ~interiorV;

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

%         sumVcol = sum(dV);
%         sumVrow = sum(dV,2);
%         sumV = sum(sumVcol);

        TdV = T*dV;
%         weightMatrixTimesVvectorized(net, dV, ...
%             sumVrow, sumVcol, sumV);

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

%         sumVcol = sum(V);
%         sumVrow = sum(V,2);
%         sumV = sum(sumVcol);
% 
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

function net = computeSolution(net,V,iter)

    V(V > 1 - power(10, -1 * net.setting.e)) = 1; %FIXME to be removed?
    V(V < power(10, -1 * net.setting.e)) = 0;

    if isequal(net.problemParameters.R*V,net.problemParameters.b)
        net.results.validSolution = true;
        net.results.exitFlag = 1;
    else
        % Error produced during simulation
        if iter > net.setting.maxIter % Max iterations reached
            net.results.exitFlag = 0;
        else
            net.results.exitFlag = -1;
        end
        net.results.validSolution = false;
    end
    net.results.x = V;
	net.results.fval = 0.5*net.results.x'*net.problemParameters.P*net.results.x + net.problemParameters.q'*net.results.x;

    net.results.itersReached = iter;

end
