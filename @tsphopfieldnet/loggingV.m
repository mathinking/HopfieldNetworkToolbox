function loggingV(iter,V,dU)

%     toXLS = cell(size(V,1) + 3,2*size(V,2) + 2);
%     toXLS(2,1) = num2cell(iter);
%     toXLS(3:end-1,2:size(V,2)+1) = num2cell(V);
%     toXLS(3:end-1,size(V,2)+3:end) = num2cell(dU);
%     xlswrite('test',toXLS,1,['A',num2str(size(toXLS,1)*(iter-1)+1)])
    
    if iter == 1
        Vlog = nan([size(V),2000]);
        dUlog = nan([size(dU),2000]);
        save('test','Vlog','dUlog','-v7.3');
    end
    data = matfile('test','Writable',true);
    data.Vlog(:,:,iter) = V;
    data.dUlog(:,:,iter) = dU;
    
end