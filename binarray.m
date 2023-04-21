% Simple function for binning array values in to bins of size n
%   Similar to `decimate`, but works on GPU
%   x = binarray(x, n, fun)
function x = binarray(x, n, fun)
    if n <= 1   % No decimation
        return;
    end
    
    N = 1;
    % Set default for 'fun' input
    if ~exist('fun', 'var')
        if isgpuarray(x)
            fun = @sum;
            N = n;
        else
            fun = @mean;
        end
    end
    
    % Create bins
    bins = reshape(repmat(1:ceil(numel(x)/n), n, 1), [], 1);
    bins = bins(1:numel(x));
    
    % Select the minimum from each coordinate
    x = accumarray(bins, x(:)/N, [ceil(numel(x)/n), 1], fun, NaN);
end

