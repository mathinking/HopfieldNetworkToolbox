function setSimFcn(net,simFcn)
   
    options = tsphopfieldnetOptions('SimFcn', simFcn);   
    net.SimFcn = options.SimFcn;
    
end

