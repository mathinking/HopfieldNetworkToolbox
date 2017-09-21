mainFolder = fileparts(fileparts(which('hopfieldnetwork')));
customXsl = fullfile(mainFolder,'help_source','mxdom2mathjax.xsl');

addpath(mainFolder)
addpath(fullfile(mainFolder,'chn'));
addpath(genpath(fullfile(mainFolder,'data','TSPFiles')));

helpFolder = fullfile(mainFolder,'help');
helpSourceFilesLocation = fullfile(mainFolder,'help_source');

addpath(helpFolder)
addpath(fullfile(helpSourceFilesLocation))

optionsPublish = struct('format','html','outputDir',fullfile(mainFolder,'help','html'),'stylesheet',customXsl);

publish(fullfile(helpSourceFilesLocation,'chn_product_page'),optionsPublish);
publish(fullfile(helpSourceFilesLocation,'chn_about'),optionsPublish);
publish(fullfile(helpSourceFilesLocation,'chn_concepts'),optionsPublish);
publish(fullfile(helpSourceFilesLocation,'chn_release_notes'),optionsPublish);
publish(fullfile(helpSourceFilesLocation,'chn_system_requirements'),optionsPublish);
publish(fullfile(helpSourceFilesLocation,'chn_users_guide'),optionsPublish);
publish(fullfile(helpSourceFilesLocation,'chn_users_guide_hopfieldnet'),optionsPublish);
publish(fullfile(helpSourceFilesLocation,'chn_users_guide_tsplib'),optionsPublish);
publish(fullfile(helpSourceFilesLocation,'chn_users_guide_howToUseTsphopfieldnetOptions'),optionsPublish);
publish(fullfile(helpSourceFilesLocation,'chn_users_guide_howToUsehopfieldnetOptions'),optionsPublish);
publish(fullfile(helpSourceFilesLocation,'chn_users_guide_TSP'),optionsPublish);
publish(fullfile(helpSourceFilesLocation,'chn_users_guide_SimulinkModels'),optionsPublish);
publish(fullfile(helpSourceFilesLocation,'chn_references'),optionsPublish);
publish(fullfile(helpFolder,'Example_tspUsingRegularPolygons'));
publish(fullfile(helpFolder,'Example_tspUsingTSPLIB'));
publish(fullfile(helpFolder,'Example_tspUsingCoords'));
publish(fullfile(helpFolder,'Example_tspUsingDistance'));
publish(fullfile(helpSourceFilesLocation,'chn_users_guide_improve'),optionsPublish);
publish(fullfile(helpFolder,'Example_tspReducingC'));
publish(fullfile(helpFolder,'Example_tspSaddlePoint'));
publish(fullfile(helpSourceFilesLocation,'chn_users_guide_improveHybrid'),optionsPublish);
publish(fullfile(helpFolder,'Example_tspDivideConquer'));
publish(fullfile(helpFolder,'Example_tsp2opt'));
publish(fullfile(helpSourceFilesLocation,'chn_TSP_APP'),optionsPublish);
publish(fullfile(helpFolder,'Example_GQKPusingCHN'));
publish(fullfile(helpFolder,'Example_SimulinkDiscrete'));
publish(fullfile(helpFolder,'Example_SimulinkContinuous'));

builddocsearchdb(fullfile(pwd,helpFolder,'html'));

rmpath(fullfile(helpSourceFilesLocation))

close all;
