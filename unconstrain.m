%% Inverse Hyperbolic Constraint Function
function x = unconstrain(xp, low, high, alpha)
    if ~exist('alpha', 'var')
        alpha = 10*ones(size(low));
    end
    if ~exist('high', 'var')
        high = 1;
    end
    if ~exist('low', 'var')
        low = 1;
    end
    
    ldim = find(size(low) == min(size(low)), 1);
    range = abs(high - low);
    offset = mean(cat(ldim, low, high), ldim);
    
    x = -atanh( 2*(offset - xp) ./ range ) .* alpha * pi;
end
