function setTrainFcn(net,trainFcn)

    options = hopfieldnetOptions('TrainFcn', trainFcn);   
    net.TrainFcn = options.TrainFcn;
    
end
