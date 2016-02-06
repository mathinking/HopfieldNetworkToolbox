function setSimFcn(net,simFcn)

    assert(strcmp('euler',simFcn) || strcmp('talavan-yanez',simFcn) || strcmp('talavan-yanez-neighbors',simFcn), ...        
        'Simulation function must be ''euler'' or ''talavan-yanez'' or ''talavan-yanez-neighbors''');
    net.simFcn = simFcn;
    
end

