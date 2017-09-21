function energyplot(net,varargin)
    if nargin > 1
        ax = varargin{1};
    else
        ax = [];
    end
    if net.Results.ValidSolution
        if isempty(ax)
            ax = gca;
        end
        toBePlotted = 1:net.Results.ItersReached-1;
        plot(ax,net.Results.Time(toBePlotted), net.Results.Energy(toBePlotted), '.b-');
        title(ax,'Energy Function');
        if nargin == 1
            xlabel(ax,'Time');
        end
        axis(ax,'tight')
        ylim([min(net.Results.Energy(toBePlotted)), max(net.Results.Energy(toBePlotted))])
    else
        error('tsphopfieldnet:energyplot','Unvalid energy results. Possibly no path has been found or training has not taken place yet.');
    end
end
