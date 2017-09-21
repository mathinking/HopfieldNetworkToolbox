function energyplot(net,varargin)
    if nargin > 1
        ax = varargin{1};
    else
        ax = [];
    end
    if net.Results.ValidPath
        if isempty(ax)
            ax = axes;
        end
        toBePlotted = ~isnan(net.Results.Energy) & net.Results.Time < 10^5;
        if ~isempty(varargin)
            toBePlotted(1) = 0; % Not plotting t = 0
        end
        plot(ax,net.Results.Time(toBePlotted), net.Results.Energy(toBePlotted), '.b-');
        title(ax,'Energy Function');
        if nargin == 1
            xlabel(ax,'Time');
        end
        axis(ax,'tight')
        ylim(ax, [min(net.Results.Energy(toBePlotted)), max(net.Results.Energy(toBePlotted))])      
    else
        error('tsphopfieldnet:energyplot','Unvalid energy results. Possibly no path has been found or training has not taken place yet.');
    end
end
