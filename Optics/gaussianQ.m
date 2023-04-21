%% Helper function to calculate Gaussian Complex Beam Parameter
% q = gaussianQ(lambda, w0[, d])
function q = gaussianQ(lambda, w0, d)
    if ~exist('d', 'var')
        d = 0;
    end
    q = d + 1i*pi*w0^2/lambda;
end
