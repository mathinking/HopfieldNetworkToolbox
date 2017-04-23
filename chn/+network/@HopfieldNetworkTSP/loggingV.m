function loggingV(checkpointFilename,maxIter,iter,N,V,dU)
   
    if iter == 1
        Vlog  = nan([N,N,maxIter]); %#ok<NASGU>
        dUlog = nan([N,N,maxIter]); %#ok<NASGU>
        save(checkpointFilename,'Vlog','dUlog','-v7.3');
        V  = zeros(N);
        dU = zeros(N);
    end
    data = matfile(checkpointFilename,'Writable',true);
    if isa(V,'gpuArray')
        V = gather(V);
        dU = gather(dU);
    end
    data.Vlog(:,:,iter)  = V;
    data.dUlog(:,:,iter) = dU;
    
end
