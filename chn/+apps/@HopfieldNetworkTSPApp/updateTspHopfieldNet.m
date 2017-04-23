function updateTspHopfieldNet(app,C)
    setTrainParam(app.net, 'C', C);
    train(app.net);
    
    trainParams = getTrainParam(app.net);

    app.parameterAEdit.String = num2str(trainParams.A);
    app.parameterBEdit.String = num2str(trainParams.B);
    app.parameterDEdit.String = num2str(1/trainParams.dUaux);
    app.parameterNpEdit.String = num2str(trainParams.Np);
end
