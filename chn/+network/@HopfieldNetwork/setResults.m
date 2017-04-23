function setResults(net,property,value)

    if isfield(net.Results,property)
        net.Results.(property) = value;
    else
        error(['There is no results property named ', property]);
    end
 
end
