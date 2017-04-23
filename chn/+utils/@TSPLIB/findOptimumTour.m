function [tourLength,optimumTour] = findOptimumTour(problem)

    fid = fopen([problem.name,'.opt.tour'],'rt');
    if fid < 0
        fprintf(['Optimum tour file not available for ',...
            'problem %s.\n'],problem.name);
        tourLength = NaN;
        optimumTour = NaN;
        return;
    end    
    positioned = false;

    while ~positioned
       str = fgetl(fid);
       if strcmp(str,'TOUR_SECTION') || feof(fid)
           positioned = true;
       end
    end

    optimumTour = textscan(fid, '%f');
    optimumTour = optimumTour{1}'; 

    fid = fclose(fid);      %#ok<NASGU>

    optimumTour(end) = optimumTour(1);
    tourLength = problem.computeTour(problem.d,optimumTour);            

end
