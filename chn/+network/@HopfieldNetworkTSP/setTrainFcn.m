function setTrainFcn(net,trainFcn)

    options = tsphopfieldnetOptions('TrainFcn', trainFcn);   
    net.TrainFcn = options.TrainFcn;
    
end
