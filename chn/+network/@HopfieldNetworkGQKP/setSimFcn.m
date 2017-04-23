function setSimFcn(net,simFcn)

    options = hopfieldnetOptions('SimFcn', simFcn);   
    net.SimFcn = options.SimFcn;

end

