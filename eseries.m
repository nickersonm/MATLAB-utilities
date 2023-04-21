% Full E series Em; preferred numbers on logarithmic scale
%   https://en.wikipedia.org/wiki/E_series_of_preferred_numbers
function s = eseries(m)
    if m <= 24; rN = 1; else; rN = 2; end
    s = round((10.^(0:(m-1))).^(1/m), rN);
end
