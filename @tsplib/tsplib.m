classdef tsplib < handle
    
    properties (GetAccess = public, SetAccess = private)
        name;
        nCities
        coords;
        type;
        d;
    end
    
    methods %(Sealed = true)       
    	function problem = tsplib(name)
            if nargin ~= 0 % Allow nargin == 0 syntax
                problem(length(name),1) = tsplib;
                for i = 1:length(name)

                    problem(i).name = name{i}; %#ok<*AGROW>

                    fid = fopen([problem(i).name,'.tsp'],'rt');

                    positioned = false;

                    while ~positioned
                        str = fgetl(fid);
                        if strncmp(str,'DIMENSION', 9)
                            problem(i).nCities = str2double(regexprep(...
                                str(length('DIMENSION')+1:end),...
                                {':',' '},''));
                        end
                        if strncmp(str, 'EDGE_WEIGHT_TYPE',16)
                            problem(i).type = regexprep(...
                                str(length('EDGE_WEIGHT_TYPE')+1:end),...
                                {':',' '},'');
                        end
                        if strcmp(problem(i).type,'EXPLICIT')
                            if strncmp(str,'EDGE_WEIGHT_FORMAT', 18)
                                matrixType = regexprep(...
                                    str(length('EDGE_WEIGHT_FORMAT')+1:end),...
                                    {':',' '},'');
                            end
                            if strncmp(str, 'EDGE_WEIGHT_SECTION',19)
                                positioned = true;
                            end    
                        else
                            if strcmp(str,'NODE_COORD_SECTION') || ...
                                strcmp(str,'DISPLAY_DATA_SECTION') || feof(fid)
                                positioned = true;
                            end
                        end
                    end

                    if ~strcmp(problem(i).type,'EXPLICIT')
                        problem(i).coords = textscan(fid, '%*u%f%f%*[^\n]');
                        problem(i).coords = [problem(i).coords{:}];
                    else
                        dataCoords = textscan(fid, '%f');
                        dataCoords = [dataCoords{:}]; 

                        if strcmp(matrixType, 'LOWER_DIAG_ROW')
                            problem(i).d = ones(problem(i).nCities);
                            problem(i).d(logical(tril(problem(i).d))') = dataCoords;
                            problem(i).d = problem(i).d-tril(problem(i).d)+triu(problem(i).d)';
                        elseif strcmp(matrixType, 'UPPER_ROW')    
                            problem(i).d = ones(problem(i).nCities);
                            problem(i).d(logical(triu(problem(i).d,1))') = dataCoords;
                            problem(i).d = problem(i).d + triu(problem(i).d',1) - ...
                                triu(ones(problem(i).nCities));
                        elseif strcmp(matrixType, 'FULL_MATRIX')
                            problem(i).d = reshape(dataCoords,...
                                problem(i).nCities,problem(i).nCities);
                        elseif strcmp(matrixType, 'UPPER_DIAG_ROW')
                            problem(i).d = ones(problem(i).nCities);
                            problem(i).d(logical(tril(problem(i).d))) = dataCoords;
                            problem(i).d = problem(i).d - triu(problem(i).d);
                            problem(i).d = problem(i).d + problem(i).d';
                        end
                    end

                    fid = fclose(fid);      %#ok<NASGU>

                    if strcmp(problem(i).type,'GEO')
                        problem(i).coords = tsplib.convert2LatLon(problem(i).coords); %#TODO Should it be in the object instantation?
                        problem(i).d = tsphopfieldnet.computeDistance(problem(i).coords,problem(i).type);
                    elseif strcmp(problem(i).type,'EUC_2D')   
                        problem(i).d = tsphopfieldnet.computeDistance(problem(i).coords,problem(i).type);
                    elseif strcmp(problem(i).type,'CEIL_2D')
                        problem(i).d = tsphopfieldnet.computeDistance(problem(i).coords,problem(i).type);
                    elseif strcmp(problem(i).type,'ATT')   
                        problem(i).d = tsphopfieldnet.computeDistance(problem(i).coords,problem(i).type);
                    elseif strcmp(problem(i).type,'EUC') 
                        problem(i).d = tsphopfieldnet.computeDistance(problem(i).coords,problem(i).type);
                    end
                end
            end
        end
    end   

    % --- Methods definitions --- %    
    methods (Hidden = false, Access = public)
        plot(problem);
        [tourLength,optimumTour] = findOptimumTour(problem);
    end
    
    methods (Static = true, Hidden = false, Access = public)
        problemNames = problemNames(hasOptimumTourFile);
    end
    
    methods (Static = true, Hidden = true, Access = public)
        dataCoords = convert2LatLon(dataCoords);
        tourLength = computeTour(d,optimumTour);
    end
    
end
