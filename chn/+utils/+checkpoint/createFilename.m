function checkpointFilename = createFilename(checkpointPath,simFcn)

    warning('off','MATLAB:MKDIR:DirectoryExists');
    mkdir(checkpointPath);
    warning('on','MATLAB:MKDIR:DirectoryExists');
    checkpointFilename = ['checkpoint_',simFcn,'_',...
        datestr(datetime(now,'ConvertFrom','datenum'),...
                'yyyy_mm_dd_HH_MM_SS'),'.mat'];

end