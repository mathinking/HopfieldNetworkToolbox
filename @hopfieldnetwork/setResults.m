function setResults(net,property,value)

    if isfield(net.results,property)
        net.results.(property) = value;
    else
        error(['There is no results property named ', property]);
    end
 
end
