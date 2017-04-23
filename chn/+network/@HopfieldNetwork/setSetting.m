function setSetting(net,property,value)

    if isfield(net.setting,property)
        hopfieldnetwork.isValidSetting(property,value)
        
        % Transfer function must be either satlin or tanh
        if strcmp(property,'transferFcn') || strcmp(property,'invTransferFcn')
            if strcmp(value,'satlin') || strcmp(value,'invsatlin') 
                net.setting.transferFcn = @(u) net.satlin(u, net.setting.u0);
                net.setting.invTransferFcn = @(v) net.invsatlin(v, net.setting.u0);
            else
                net.setting.transferFcn = @(u) 0.5 * (1 + tanh(u ./ options.setting.u0));
                net.setting.invTransferFcn = @(v) (net.setting.u0 / 2) * log(v./(1-v));
            end
        else
        	net.setting.(property) = value;
        end
    else
        error('hopfieldNetwork:unvalidSetting',['There is no setting named ', property]);
    end
 
end
