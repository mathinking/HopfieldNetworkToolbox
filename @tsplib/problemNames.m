function problems = problemNames(hasOptimumTourFile, saveMATfile)

filedir = fileparts(which('setup_hopfieldNetwork.m'));
TSPFilesDir = 'TSPFiles';
TSPToursDir = 'TSPTours';

if hasOptimumTourFile
    OptimTourTSPFiles = cellstr(ls(fullfile(filedir,TSPFilesDir,TSPToursDir,'*.tour')));

    problemsTour = regexprep(OptimTourTSPFiles,'.opt.tour','');
    TSPFiles = cellstr(ls(fullfile(filedir,TSPFilesDir,'*.tsp')));
    problems = regexprep(TSPFiles,'.tsp','');
    problems = intersect(problems,problemsTour);
    
else
    TSPFiles = cellstr(ls(fullfile(filedir,TSPFilesDir,'*.tsp')));
    problems = regexprep(TSPFiles,'.tsp','');
end

if saveMATfile
    save(fullfile(filedir,TSPFilesDir,'TSPLIBproblems.mat'),'problems')
end