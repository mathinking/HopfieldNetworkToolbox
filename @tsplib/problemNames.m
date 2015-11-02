function problems = problemNames(hasOptimumTourFile)

TSPFilesDir = 'TSPFiles';
TSPToursDir = 'TSPTours';

if hasOptimumTourFile
    % TSPFiles = cellstr(ls(fullfile(pwd,TSPFilesDir,'*.tsp')));
    OptimTourTSPFiles = cellstr(ls(fullfile(pwd,TSPFilesDir,TSPToursDir,'*.tour')));

    problems = regexprep(OptimTourTSPFiles,'.opt.tour','');
else
    TSPFiles = cellstr(ls(fullfile(pwd,TSPFilesDir,'*.tsp')));
    problems = regexprep(TSPFiles,'.tsp','');
end

% save(fullfile(pwd,TSPFilesDir,'TSPLIBproblems.mat'),'problemNames')
