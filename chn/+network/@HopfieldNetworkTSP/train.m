function net = train(net)
        
    if isfield(net.TrainParam,'dUaux') && net.TrainParam.dU == 1 %Prepare for retraining
        net.Cities.DistanceMatrix = net.TrainParam.dUaux * net.Cities.DistanceMatrix;
    end
    
    if strcmp(net.TrainFcn,'trainty')
        net.TrainParam.dU = max(net.Cities.DistanceMatrix(triu(net.Cities.DistanceMatrix~= 0))); %max(net.Cities.DistanceMatrix(triu(true(net.TrainParam.N),1)));
        
        if strcmp(net.SimFcn,'talavan-yanez')
            net.Cities.DistanceMatrix = net.Cities.DistanceMatrix / net.TrainParam.dU;
            net.TrainParam.dUaux = net.TrainParam.dU;
            net.TrainParam.dU = 1;
        end
        net.TrainParam.dL = min(net.Cities.DistanceMatrix(triu(net.Cities.DistanceMatrix~= 0))); %min(net.Cities.DistanceMatrix(triu(true(net.TrainParam.N),1)));

        net.TrainParam.D = 1/net.TrainParam.dU;
        net.TrainParam.rho = net.TrainParam.dL/net.TrainParam.dU;

        net.TrainParam.Np = net.TrainParam.N - net.TrainParam.K + 3 / net.TrainParam.C;

        if net.TrainParam.K == 0
            net.TrainParam.B = 3 + net.TrainParam.C;
            net.TrainParam.A = net.TrainParam.B - net.TrainParam.rho;
        else
            net.TrainParam.A = 3 + net.TrainParam.C;
            net.TrainParam.B = net.TrainParam.A + net.TrainParam.rho;
        end        

        net.TrainParam = orderfields(net.TrainParam);
    else
        net.TrainFcn = 'trainty';
        error('tsphopfieldnet:unvalidTrainFcn', ['Not Valid Training Function. Setting to ''trainty''. ',...
            'Please, train the network using train(net)']);
    end
    
end
