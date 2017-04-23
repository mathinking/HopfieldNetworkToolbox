function net = train(net)
        
    if isfield(net.trainParam,'dUaux') && net.trainParam.dU == 1 %Prepare for retraining
        net.cities.d = net.trainParam.dUaux * net.cities.d;
    end
    
    if strcmp(net.trainFcn,'trainty')
        net.trainParam.dU = max(net.cities.d(triu(net.cities.d~= 0))); %max(net.cities.d(triu(true(net.trainParam.N),1)));
        
        if strcmp(net.simFcn,'talavan-yanez')
            net.cities.d = net.cities.d / net.trainParam.dU;
            net.trainParam.dUaux = net.trainParam.dU;
            net.trainParam.dU = 1;
        end
        net.trainParam.dL = min(net.cities.d(triu(net.cities.d~= 0))); %min(net.cities.d(triu(true(net.trainParam.N),1)));

        net.trainParam.D = 1/net.trainParam.dU;
        net.trainParam.rho = net.trainParam.dL/net.trainParam.dU;

        net.trainParam.Np = net.trainParam.N - net.trainParam.K + 3 / net.trainParam.C;

        if net.trainParam.K == 0
            net.trainParam.B = 3 + net.trainParam.C;
            net.trainParam.A = net.trainParam.B - net.trainParam.rho;
        else
            net.trainParam.A = 3 + net.trainParam.C;
            net.trainParam.B = net.trainParam.A + net.trainParam.rho;
        end        

        net.trainParam = orderfields(net.trainParam);
    else
        net.trainFcn = 'trainty';
        error('tsphopfieldnet:unvalidTrainFcn', ['Not Valid Training Function. Setting to ''trainty''. ',...
            'Please, train the network using train(net)']);
    end
    
end
