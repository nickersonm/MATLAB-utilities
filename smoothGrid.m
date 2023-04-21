%% Helper function to simplify gridding
% Michael Nickerson 2019
% [rowg, colg, zg, binC] = smoothGrid(row, col, z, N)
% Turns arbitrary scattered data into a smoothed NxN grid
%   Uses smoothn(<z>, 1) for smoothing
function [rowg, colg, zg, binC] = smoothGrid(row, col, z, N, fun)
    if ~exist('fun', 'var'); fun = @min; end
    if numel(N) == 1; N = [N N]; end
    
    % Generate a grid and get the locations
    [binC, rowg, colg, rowi, coli] = histcounts2(row, col, N-1);
    
    % Select the minimum from each coordinate
    zg = accumarray([rowi coli], z, N, fun, NaN);
    
    % Smooth this out
    zg = smoothn(zg, 1);
end
