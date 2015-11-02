function energyplot(net)
    if net.results.validPath
        figure;
        toBePlotted = ~isnan(net.results.energy) & net.results.time < 10^100;
        semilogy(net.results.time(toBePlotted), net.results.energy(toBePlotted), '.b-');
        title('Energy Function');
        xlabel('Time');
        axis tight;
    else
        error('tsphopfieldnet:energyplot','Unvalid energy results. Possibly no path has been found or training has not taken place yet.');
    end
end
