function myText = cityTextGeneration(N)
% CITYTEXTGENERATION creates a non-repeated sequence of N city names
%
%   MYTEXT = CITYTEXTGENERATION(N) generates a sequence of N city names
%   using characters from A-Z. If more characters are needed, they are
%   generated using combinations of n things taken k at a time.
%   
%   See also NCHOOSEK

    enough = false;
    i = 1;
    j = 0;
    myPrevText = char(65:90)';
    myText = [];
    while ~enough
        if N > j + nchoosek(length(myPrevText),i)
            j = j + nchoosek(length(myPrevText),i);
            myText = char(myPrevText,...
                nchoosek('ABCDEFGHIJKLMNOPQRSTUVWXYZ',i+1));
            i = i+1;
        else 
            enough = true;
        end
    end
    if isempty(myText)
        myText = myPrevText(1:N,:);
    else
        myText = myText(1:N,:);
    end
    myText = cellstr(myText)';
end
