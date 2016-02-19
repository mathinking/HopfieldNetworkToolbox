function timeElapsed(app,~,event)

if strcmp(event.Type,'StartFcn')
    app.elapsedTimeSecs = tic;
else
    elapsedT = toc(app.elapsedTimeSecs);
    app.elapsedTime.String = ['Elapsed time: ', elapsedTimeDisplay(elapsedT)];
    drawnow
end

function str = elapsedTimeDisplay(dur)
% Creates elapsed time string in the format of hr/min/s

if dur >= 3600
    str = sprintf('%d hr %d min %0.1f s', floor(dur/3600), ...
        floor(mod(dur, 3600)/60), mod(dur, 60));
elseif dur >= 60
    str = sprintf('%d min %0.1f s', floor(dur/60), mod(dur, 60));
else
    str = sprintf('%0.1f s', dur);
end
