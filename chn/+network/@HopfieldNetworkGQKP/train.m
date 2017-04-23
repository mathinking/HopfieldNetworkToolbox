function net = train(net)
    net.TrainParam.alpha = 0;
    net.TrainParam.beta = zeros(net.ProblemParameters.m,1);
    net.TrainParam.eps = 0;
    net.TrainParam.gamma = zeros(net.ProblemParameters.n + net.ProblemParameters.m1);
    net.TrainParam.phi = zeros(net.ProblemParameters.m);
    
    e_k = {net.ProblemParameters.R,net.ProblemParameters.b};
    Wk0 = cell(net.ProblemParameters.m1,1);
    Wk1 = cell(net.ProblemParameters.m,1);
    Wk2 = cell(net.ProblemParameters.m,1);
    n = net.ProblemParameters.n;
    
%     for k = 1:net.ProblemParameters.m
%         if k <= net.ProblemParameters.m1
%             if k > 1
%                                 
%             end
%             Wk0{k}{1} = {e_k{1}(k,:),e_k{2}(k)}; % ei[k] > b[k]
%             Wk0{k}{2} = {e_k{1}(k,:) - R(k,n+k)}
%         else
%             
%         end
%     end
	
    % Gradient of E
%     E = cell(1,net.ProblemParameters.n+net.ProblemParameters.m1);
%     
%     for i = 1:length(E)
%         if i <= net.ProblemParameters.n
%             E{i} = @(alpha,Phi,beta,gamma)[alpha*P(i,:) + alpha*q(i)]
%         end
%     end
%     
    % Order of parameters: [alpha, Phi's,     
    net.TrainParam.alpha = 1;
    phi1 = 5;
    alpha = net.TrainParam.alpha;
 	net.TrainParam.eps = 3*alpha/2 + phi1/2;
    net.TrainParam.Phi = phi1;
    net.TrainParam.beta = -alpha/2 - phi1;
    net.TrainParam.gamma(1,1) = (2*alpha+phi1/2);
    net.TrainParam.gamma(2,2) = phi1/2 - alpha;
    
    % T  = -(alpha * P + R'*Phi*R - 2*diag(gamma));
    net.TrainParam.T  = -(net.TrainParam.alpha * net.ProblemParameters.P + net.ProblemParameters.R'*net.TrainParam.Phi*net.ProblemParameters.R - 2*net.TrainParam.gamma);
    
    % ib = -(alpha * q + R'*beta + gamma);
    net.TrainParam.ib = -(net.TrainParam.alpha * net.ProblemParameters.q + net.ProblemParameters.R'*net.TrainParam.beta + diag(net.TrainParam.gamma));
end

