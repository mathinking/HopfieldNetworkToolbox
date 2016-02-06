function setSimFcn(net,simFcn)

    assert(strcmp('euler',simFcn) || strcmp('talavan-yanez',simFcn) || strcmp('divide-conquer',simFcn), ...        
        'Simulation function must be ''euler'' or ''talavan-yanez'' or ''divide-conquer''');
    net.simFcn = simFcn;
    
end

