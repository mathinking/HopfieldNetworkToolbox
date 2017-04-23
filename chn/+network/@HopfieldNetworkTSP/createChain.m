function chain = createChain(cities)

    chain = [cities',cellstr(repmat('-',length(cities),1))]';
    chain = [chain{:}];
    chain(end) = [];
    
end
