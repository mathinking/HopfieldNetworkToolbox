function dataCoords = convert2LatLon(dataCoords)

    PI = 3.141592;
    degrees = fix(dataCoords);

    minutes = dataCoords - degrees;
    dataCoords = PI * (degrees + 5*minutes/3) / 180;

end
