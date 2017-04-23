function energyplot(net,varargin)
    if nargin > 1
        ax = varargin{1};
    else
        ax = [];
    end
    if net.results.validPath
        if isempty(ax)
            figure;
        end
        toBePlotted = ~isnan(net.results.energy) & net.results.time < 10^100;
        if isempty(ax)
            semilogy(net.results.time(toBePlotted), net.results.energy(toBePlotted), '.b-');
            title('Energy Function');
            xlabel('Time');
            axis tight;
        else
            semilogy(ax,net.results.time(toBePlotted), net.results.energy(toBePlotted), '.b-');
            title(ax,'Energy Function');
            % xlabel(ax,'Time');
            axis(ax,'tight')
        end       
    else
        error('tsphopfieldnet:energyplot','Unvalid energy results. Possibly no path has been found or training has not taken place yet.');
    end
end
