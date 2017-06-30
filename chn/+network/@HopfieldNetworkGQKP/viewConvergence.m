function fV = viewConvergence(iter,V,net,fV)

    if iter == 1
        fV = figure;
        imagesc(V,[0,1])
        axis square;
        set(get(fV,'Children'),'YTick',1:net.TrainParam.N)
        set(get(fV,'Children'),'XAxisLocation','top')
        set(get(fV,'Children'),'XTick',1:net.TrainParam.N)
        set(get(fV,'Children'),'TickDir','Out')
        set(get(get(fV,'Children'),'Children'),'tag','image')
        hold on;
        title(['iter = ',num2str(iter)]);
        for x = 1:net.TrainParam.N-1
            plot([0,net.TrainParam.N+0.5],[x+0.5,x+0.5],'Color','w','MarkerSize',2,'LineStyle','-')
        end
        for i = 1:net.TrainParam.N-1
            plot([i+0.5,i+0.5],[0,net.TrainParam.N+0.5],'Color','w','MarkerSize',2,'LineStyle','-')
        end

    else
        figure(fV);
        hChild = get(get(fV,'Children'),'Children');
        if strcmp(net.Setting.ExecutionEnvironment,'GPU')
            set(hChild(strcmp(get(hChild,'tag'),'image')),'CData',gather(V))
        else
            set(hChild(strcmp(get(hChild,'tag'),'image')),'CData',V)
        end
        set(get(get(fV,'Children'),'Title'),'String',['iter = ',num2str(iter)]);
        drawnow;
        pause(1-net.Setting.SimulationPlotPauseTime);
    end
    
end