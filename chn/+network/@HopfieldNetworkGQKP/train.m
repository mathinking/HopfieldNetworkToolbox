function net = train(net,phi,alpha,beta,eps,gamma)

    if nargin < 4
        if nargin == 3
            net.TrainParam.alpha = alpha;
        else
            net.TrainParam.alpha = 1;
        end
        if nargin > 1
            net.TrainParam.Phi = phi;
        else
            net.TrainParam.Phi = 3;        
        end

        net.TrainParam.eps = 3*net.TrainParam.alpha/2 + net.TrainParam.Phi/2;

        net.TrainParam.beta = -net.TrainParam.alpha/2 - net.TrainParam.Phi;
        net.TrainParam.gamma(1,1) = (2*net.TrainParam.alpha+net.TrainParam.Phi/2);
        net.TrainParam.gamma(2,2) = net.TrainParam.Phi/2 - net.TrainParam.alpha; 

    elseif nargin == 6
        net.TrainParam.alpha = alpha;
        net.TrainParam.Phi = phi;
        net.TrainParam.beta = beta;
        net.TrainParam.eps = eps;
        net.TrainParam.gamma = gamma;

    else
       error('HopfieldNetworkGQKP:train:MissingParameters', ...
           'Missing parameters required to compute weight matrix and biases');
    end
    
    % T  = -(alpha * P + R'*Phi*R - 2*diag(gamma));
    net.TrainParam.T  = -(net.TrainParam.alpha * net.ProblemParameters.P + net.ProblemParameters.R'*net.TrainParam.Phi*net.ProblemParameters.R - 2*net.TrainParam.gamma);
    
    % ib = -(alpha * q + R'*beta + diag(gamma));
    net.TrainParam.ib = -(net.TrainParam.alpha * net.ProblemParameters.q + net.ProblemParameters.R'*net.TrainParam.beta + diag(net.TrainParam.gamma));
end

