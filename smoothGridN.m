%% Helper function to simplify gridding
% Michael Nickerson 2019
% [[xg, yg, ...], fg] = smoothGridN([x, y, ..., f], N)
% Turns scattered data of arbitrary dimension into a smoothed grid
%   Uses smoothn(<z>, 1, 'robust') for smoothing
function [outC, outD] = smoothGridN(inCD, N)
    % Remove duplicates
    inCD = unique(inCD, 'rows');
    
    % Initialize variables
    d = size(inCD,2)-1;
    idx = NaN(size(inCD(:,1:d)));  outC = NaN(N,d);
    
    % Bin each coordinate dimension
    for i=1:d
        [idx(:,i), edgeD] = discretize(inCD(:,i), N);
        outC(:,i) = edge2center(edgeD); % Want centers
    end
    
    % Bin values
    outD = accumarray(idx, inCD(:,d), repmat(N, [1 d]), @mean, NaN);
    
    % Smooth
    outD = smoothn(outD, 2, 'robust');
    
    
    % Helper function
    function x = edge2center(x)
        x = diff(x)/2 + x(1:end-1);
    end
end
