% Wrap a title
% Michael Nickerson 2022
function tit = titlewrap(tit, len, split)
    if ~exist('split', 'var'); split = ","; end
    if strlength(tit) > len
        tParts = strsplit(tit, split);
        tit = tParts(1); l = strlength(tit);
        for tPart = tParts(2:end)
            if l + 2 + strlength(tPart) > len
                tit = tit + split + newline + tPart;
                l = strlength(tPart);
            else
                tit = tit + split + " " + tPart;
                l = l + 2 + strlength(tPart);
            end
        end
    end
end
