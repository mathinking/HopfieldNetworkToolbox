function net = addDefaultOptionValues(net, options)
   
    if isempty(options.trainFcn)
        net.trainFcn = 'traingty';
    end
    
    if isempty(options.simFcn)
        net.simFcn = 'talavan-yanez';
    end
        
end
