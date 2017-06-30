function trimToSimulatedData(checkpointPath, checkpointFileName, iter)

    data = matfile(fullfile(checkpointPath,...
        checkpointFileName),'Writable',true);
    data.Vlog(:,:,iter+1:end) = [];
    data.dUlog(:,:,iter+1:end) = [];
    
end