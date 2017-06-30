function setScheme(net,scheme)
   
    options = tsphopfieldnetOptions('Scheme', scheme);   
    net.Scheme = options.Scheme;
    
end

