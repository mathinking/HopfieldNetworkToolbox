function tourLength = computeTour(d,optimumTour)

    N = size(d,1);
    tourLength = 0;
    city = 1;
    while city <= N
        tourLength = tourLength + d(optimumTour(city),...
            optimumTour(city+1));
        city = city + 1;
    end

end
