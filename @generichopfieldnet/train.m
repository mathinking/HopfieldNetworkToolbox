function net = train(net)
    net.trainParam.alpha = 0;
    net.trainParam.beta = zeros(net.problemParameters.m,1);
    net.trainParam.eps = 0;
    net.trainParam.gamma = zeros(net.problemParameters.n + net.problemParameters.m1);
    net.trainParam.phi = zeros(net.problemParameters.m);
    
    e_k = {net.problemParameters.R,net.problemParameters.b};
    Wk0 = cell(net.problemParameters.m1,1);
    Wk1 = cell(net.problemParameters.m,1);
    Wk2 = cell(net.problemParameters.m,1);
    n = net.problemParameters.n;
    
%     for k = 1:net.problemParameters.m
%         if k <= net.problemParameters.m1
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
%     E = cell(1,net.problemParameters.n+net.problemParameters.m1);
%     
%     for i = 1:length(E)
%         if i <= net.problemParameters.n
%             E{i} = @(alpha,Phi,beta,gamma)[alpha*P(i,:) + alpha*q(i)]
%         end
%     end
%     
    % Order of parameters: [alpha, Phi's,     
    net.trainParam.alpha = 0.2;
    alpha = net.trainParam.alpha;
% 	net.trainParam.eps = 
    phi1 = 0.1;
    net.trainParam.Phi = [phi1,0;0,0];
    net.trainParam.beta = [0;0];
    net.trainParam.gamma(1,1) = 0.5*(alpha+phi1);
    net.trainParam.gamma(2,2) = 3/2*alpha + 0.5*phi1;
    
    % T  = -(alpha * P + R'*Phi*R - 2*diag(gamma));
    net.trainParam.T  = -(net.trainParam.alpha * net.problemParameters.P + net.problemParameters.R'*net.trainParam.Phi*net.problemParameters.R - 2*net.trainParam.gamma);
    
    % ib = -(alpha * q + R'*beta + gamma);
    net.trainParam.ib = -(net.trainParam.alpha * net.problemParameters.q + net.problemParameters.R'*net.trainParam.beta + diag(net.trainParam.gamma));
end

