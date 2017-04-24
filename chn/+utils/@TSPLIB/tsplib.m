classdef TSPLIB < handle
    %tsplib Define a TSPLIB problem
    %
    %   This class allows to define a TSPLIB problem.
    %
    %   tsplib properties:
    %       Name                       - Name of the TSPLIB problem.
    %       NumberOfCities             - Number of cities.
    %       Coordinates                - TSP problem coordinates
    %       DistanceType               - Type of distance defined for the 
    %                                    TSP problem.
    %       DistanceMatrix             - Distance matrix for the TSP
    %                                    problem.
    %
    %   Example:
    %       Create a TSPLIB class for the 'berlin52' TSPLIB problem:
    % 
    %       problem = tsplib({'berlin52'});
    %
    %   Example:
    %       Obtain a list of TSPLIB problems with a solution:
    % 
    %       problems = utils.TSPLIB.problemNames(true,false);
    %
    %   See also tsphopfieldnetOptions, tsphopfieldnet.
    
    properties (GetAccess = public, SetAccess = private)
        Name;
        NumberOfCities
        Coordinates;
        DistanceType;
        DistanceMatrix;
    end
    
    methods %(Sealed = true)       
    	function problem = TSPLIB(name)
            if nargin ~= 0 
                problem(length(name),1) = utils.TSPLIB;
                for i = 1:length(name)

                    problem(i).Name = name{i}; %#ok<*AGROW>

                    fid = fopen([problem(i).Name,'.tsp'],'rt');

                    positioned = false;

                    while ~positioned
                        str = fgetl(fid);
                        if strncmp(str,'DIMENSION', 9)
                            problem(i).NumberOfCities = str2double(regexprep(...
                                str(length('DIMENSION')+1:end),...
                                {':',' '},''));
                        end
                        if strncmp(str, 'EDGE_WEIGHT_TYPE',16)
                            problem(i).DistanceType = regexprep(...
                                str(length('EDGE_WEIGHT_TYPE')+1:end),...
                                {':',' '},'');
                        end
                        if strcmp(problem(i).DistanceType,'EXPLICIT')
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

                    if ~strcmp(problem(i).DistanceType,'EXPLICIT')
                        problem(i).Coordinates = textscan(fid, '%*u%f%f%*[^\n]');
                        problem(i).Coordinates = [problem(i).Coordinates{:}];
                    else
                        dataCoords = textscan(fid, '%f');
                        dataCoords = [dataCoords{:}]; 

                        if strcmp(matrixType, 'LOWER_DIAG_ROW')
                            problem(i).DistanceMatrix = ones(problem(i).NumberOfCities);
                            problem(i).DistanceMatrix(logical(tril(problem(i).DistanceType))') = dataCoords;
                            problem(i).DistanceMatrix = problem(i).DistanceMatrix-tril(problem(i).DistanceMatrix)+triu(problem(i).DistanceMatrix)';
                        elseif strcmp(matrixType, 'UPPER_ROW')    
                            problem(i).DistanceMatrix = ones(problem(i).NumberOfCities);
                            problem(i).DistanceMatrix(logical(triu(problem(i).DistanceMatrix,1))') = dataCoords;
                            problem(i).DistanceMatrix = problem(i).DistanceMatrix + triu(problem(i).DistanceMatrix',1) - ...
                                triu(ones(problem(i).NumberOfCities));
                        elseif strcmp(matrixType, 'FULL_MATRIX')
                            problem(i).DistanceMatrix = reshape(dataCoords,...
                                problem(i).NumberOfCities,problem(i).NumberOfCities);
                        elseif strcmp(matrixType, 'UPPER_DIAG_ROW')
                            problem(i).DistanceMatrix = ones(problem(i).NumberOfCities);
                            problem(i).DistanceMatrix(logical(tril(problem(i).DistanceMatrix))) = dataCoords;
                            problem(i).DistanceMatrix = problem(i).DistanceMatrix - triu(problem(i).DistanceMatrix);
                            problem(i).DistanceMatrix = problem(i).DistanceMatrix + problem(i).DistanceMatrix';
                        end
                    end

                    fid = fclose(fid);      %#ok<NASGU>

                    if strcmp(problem(i).DistanceType,'GEO')
                        problem(i).Coordinates = utils.TSPLIB.convert2LatLon(problem(i).Coordinates); %#TODO Should it be in the object instantation?
                        problem(i).DistanceMatrix = network.HopfieldNetworkTSP.computeDistance(problem(i).Coordinates,problem(i).DistanceType);
                    elseif strcmp(problem(i).DistanceType,'EUC_2D')   
                        problem(i).DistanceMatrix = network.HopfieldNetworkTSP.computeDistance(problem(i).Coordinates,problem(i).DistanceType);
                    elseif strcmp(problem(i).DistanceType,'CEIL_2D')
                        problem(i).DistanceMatrix = network.HopfieldNetworkTSP.computeDistance(problem(i).Coordinates,problem(i).DistanceType);
                    elseif strcmp(problem(i).DistanceType,'ATT')   
                        problem(i).DistanceMatrix = network.HopfieldNetworkTSP.computeDistance(problem(i).Coordinates,problem(i).DistanceType);
                    elseif strcmp(problem(i).DistanceType,'EUC') 
                        problem(i).DistanceMatrix = network.HopfieldNetworkTSP.computeDistance(problem(i).Coordinates,problem(i).DistanceType);
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
        problemNames = problemNames(hasOptimumTourFile, saveMATfile);
    end
    
    methods (Static = true, Hidden = true, Access = public)
        dataCoords = convert2LatLon(dataCoords);
        tourLength = computeTour(d,optimumTour);
    end
    
end
