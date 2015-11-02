function [modifiedDistance,Ng] = neighbourDistance(net, tau_or_p)

if isscalar(tau_or_p)
    if tau_or_p > 1

        tau = tau_or_p;
        Ng = false(size(net.cities.d));
        for x = 1:size(net.cities.d,1)
            neighbourCities = findClosestCities(net.cities.d(x,:),tau);
            Ng(x,neighbourCities) = true;
            Ng(neighbourCities,x) = true;
        end

        modifiedDistance = zeros(size(Ng));
%         modifiedDistance(Ng == 0 & ~logical(eye(size(net.cities.d)))) = net.trainParam.dU; %net.cities.d(Ng == 0 & ~logical(eye(size(net.cities.d))));%Dmax;
%         modifiedDistance(Ng == 0 & ~logical(eye(size(net.cities.d)))) = net.trainParam.dU - net.trainParam.rho./net.cities.d(Ng == 0 & ~logical(eye(size(net.cities.d))));
        modifiedDistance(Ng == 0 & ~logical(eye(size(net.cities.d)))) = net.trainParam.dU - (net.trainParam.dL./net.cities.d(Ng == 0 & ~logical(eye(size(net.cities.d))))-net.trainParam.rho);
        modifiedDistance(Ng == 1) = net.cities.d(Ng == 1);%net.trainParam.dL;%net.cities.d(Ng == 1);
%         modifiedDistance(Ng == 1) = net.trainParam.dL + (net.trainParam.dL./net.cities.d(Ng == 1 & ~logical(eye(size(net.cities.d))))-net.trainParam.rho);
%         modifiedDistance(Ng == 0) = net.cities.d(Ng == 0);
    else
        p = tau_or_p;
        alpha = p * net.trainParam.dU + (1-p) * net.trainParam.dL;
        modifiedDistance = net.cities.d;
        modifiedDistance( net.cities.d > alpha & ~logical(eye(net.trainParam.N)) ) = net.trainParam.dU;

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

    alpha = p * net.trainParam.dU + (1-p) * net.trainParam.dL;
    modifiedDistance = zeros(net.trainParam.N);
	modifiedDistance(net.cities.d <= alpha) = net.cities.d(net.cities.d <= alpha);
    
    Ng = net.cities.d <= alpha;
    
    for x = 1:size(net.cities.d,1)
        if all(modifiedDistance(x,:) == 0)
            neighbourCities = findClosestCities(net.cities.d(x,:),tau);
            Ng(x,neighbourCities) = true;
            Ng(neighbourCities,x) = true;
        end
    end
    modifiedDistance(Ng == 1 & net.cities.d > alpha) = net.cities.d(Ng == 1 & net.cities.d > alpha); 
    modifiedDistance(Ng == 0) = net.trainParam.dU;
    
end

end

function cities = findClosestCities(distances,nCities)

    [~,pos] = sort(distances);
    pos(1) = [];
    cities = pos(1:nCities);

end
