function setTrainFcn(net,trainFcn)

	assert(strcmp('trainty',trainFcn), ...
        'Training function (parametrization) must be ''trainty''');
    net.trainFcn = trainFcn;
    
end
