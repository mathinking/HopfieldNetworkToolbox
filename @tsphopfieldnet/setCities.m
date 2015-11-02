function setCities(net,property,value)

    if isfield(net.cities,property)
        tsphopfieldnet.isValidSetting(property,value)
        
        net.cities.(property) = value;
    else
        error('hopfieldNetwork:unvalidSetting',['There is no cities property named ', property]);
    end

end