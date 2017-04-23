function setSetting(net,property,value)

    if isfield(net.Setting, property) && ~strcmp(property, 'InvTransferFcn')
        if isa(net, 'tsphopfieldnet')
            options = tsphopfieldnetOptions(property,value);
        else
            options = hopfieldnetOptions(property,value);
        end
        
        if strcmp(property, 'TransferFcn') % Modifying TransferFcn automatically changes InvTransferFcn
            if strcmp(value, 'satlin')
                net.Setting.(property)     = @(u) net.satlin(u, net.Setting.U0);
                net.Setting.InvTransferFcn = @(v) net.invsatlin(v, net.Setting.U0);
            else % tanh
                net.Setting.(property) = @(u) 0.5 * (1 + tanh(u ./ net.Setting.U0));
                net.Setting.InvTransferFcn = @(v) (net.Setting.U0 / 2) * log(v./(1-v));                
            end
        else
            net.Setting.(property) = options.Setting.(property);
        end
        
    elseif strcmp(property, 'InvTransferFcn') 
        error('hopfieldNetwork:unvalidSetting',['', property,''' can only be modified by changing ''TransferFcn''']);
    else
        error('hopfieldNetwork:unvalidSetting',['There is no Setting named ', property]);
    end

end