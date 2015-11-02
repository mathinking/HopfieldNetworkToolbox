% This script performs the files setup for the TSP Problem. The data was
% downloaded from TSPLIB
% http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/

TSPFilesDir = 'TSPFiles';
TSPToursDir = 'TSPTours';
untar('ALL_tsp.tar.gz',TSPFilesDir)
gunzip([TSPFilesDir,'/*.gz'],TSPFilesDir);
delete([TSPFilesDir,'/*.gz']);

mkdir(TSPFilesDir,TSPToursDir)
movefile([TSPFilesDir,'\*.tour'],[TSPFilesDir,'\',TSPToursDir])
