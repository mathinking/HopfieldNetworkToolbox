function problems = problemNames(hasOptimumTourFile, saveMATfile)

filedir = fileparts(fileparts(which('tsphopfieldnet.m')));
TSPFilesDir = 'TSPFiles';
TSPToursDir = 'TSPTours';

if hasOptimumTourFile
    OptimTourTSPFiles = dir(fullfile(filedir,'data',TSPFilesDir,TSPToursDir,'*.tour'));
    OptimTourTSPFiles = {OptimTourTSPFiles.name}';
    
    problemsTour = regexprep(OptimTourTSPFiles,'.opt.tour','');
    TSPFiles = dir(fullfile(filedir,'data',TSPFilesDir,'*.tsp'));
    TSPFiles = {TSPFiles.name}';
    
    problems = regexprep(TSPFiles,'.tsp','');
    problems = intersect(problems,problemsTour);
    
else
    TSPFiles = dir(fullfile(filedir,'data',TSPFilesDir,'*.tsp'));
    TSPFiles = {TSPFiles.name}';
    problems = regexprep(TSPFiles,'.tsp','');
end

if saveMATfile
    save(fullfile(filedir,'data',TSPFilesDir,'TSPLIBproblems.mat'),'problems')
end
