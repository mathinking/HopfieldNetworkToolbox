function [modifiedDistance,Ng] = neighbourDistance(net, tau_or_p)

if isscalar(tau_or_p)
    if tau_or_p > 1

        tau = tau_or_p;
        Ng = false(size(net.Cities.DistanceMatrix));
        for x = 1:size(net.Cities.DistanceMatrix,1)
            neighbourCities = findClosestCities(net.Cities.DistanceMatrix(x,:),tau);
            Ng(x,neighbourCities) = true;
            Ng(neighbourCities,x) = true;
        end

        modifiedDistance = zeros(size(Ng));
        modifiedDistance(Ng == 0 & ~logical(eye(size(net.Cities.DistanceMatrix)))) = ...
            net.TrainParam.dU - (net.TrainParam.dL./net.Cities.DistanceMatrix(Ng == 0 & ~logical(eye(size(net.Cities.DistanceMatrix))))-net.TrainParam.rho);
        modifiedDistance(Ng == 1) = net.Cities.DistanceMatrix(Ng == 1);
    else
        p = tau_or_p;
        alpha = p * net.TrainParam.dU + (1-p) * net.TrainParam.dL;
        modifiedDistance = net.Cities.DistanceMatrix;
        modifiedDistance( net.Cities.DistanceMatrix > alpha & ~logical(eye(net.TrainParam.N)) ) = net.TrainParam.dU;

        Ng = modifiedDistance <= alpha;

    end
else
    if ~any(tau_or_p < 1) || ~any(tau_or_p >= 1) || numel(tau_or_p) ~= 2
        error(['Please check the values given to tau. It must be a 1 or 2 element vector ',...
            'containing the number of neighbours for points with large entropy ',...
            'and distance radius for points with low entropy'])
    end
    tau = tau_or_p(tau_or_p >= 1);
    p   = tau_or_p(tau_or_p <  1);

    alpha = p * net.TrainParam.dU + (1-p) * net.TrainParam.dL;
    modifiedDistance = zeros(net.TrainParam.N);
	modifiedDistance(net.Cities.DistanceMatrix <= alpha) = net.Cities.DistanceMatrix(net.Cities.DistanceMatrix <= alpha);
    
    Ng = net.Cities.DistanceMatrix <= alpha;
    
    for x = 1:size(net.Cities.DistanceMatrix,1)
        if all(modifiedDistance(x,:) == 0)
            neighbourCities = findClosestCities(net.Cities.DistanceMatrix(x,:),tau);
            Ng(x,neighbourCities) = true;
            Ng(neighbourCities,x) = true;
        end
    end
    modifiedDistance(Ng == 1 & net.Cities.DistanceMatrix > alpha) = net.Cities.DistanceMatrix(Ng == 1 & net.Cities.DistanceMatrix > alpha); 
    modifiedDistance(Ng == 0) = net.TrainParam.dU;
    
end

end

function cities = findClosestCities(distances,nCities)

    [~,pos] = sort(distances);
    pos(1) = [];
    cities = pos(1:nCities);

end
