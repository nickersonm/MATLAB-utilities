%% Compute convolution via FFT
% Michael Nickerson 2019
% Replicates behavior of 'convn' via FFT for faster performance
% Accepts shapes 'same' and 'full' (default) only
function w = fftconv(u, v, shape)
    if ~exist('shape', 'var')
        shape = 'full';
    end
    
    n = 2*max(size(u), size(v)) - 1;
    w = ifftn( fftn(u, n) .* fftn(v, n) );
    
    if strcmpi(shape, 'same')
        % TODO: Impove this logic
        inds = cell(size(size(w)));
        for i = 1:numel(size(w))
            inds{i} = floor( size(u,i)/2+1 ):floor( 3*size(u,i)/2 );
        end
        
        w = w(inds{:});
    end
end
%% Tests
% i1=rand(100,1); i2=rand(100,1); max(abs(1-conv(i1, i2)./fftconv(i1,i2)))
% i1=rand(100); i2=rand(100); max(max(abs(1-convn(i1, i2, 'same')./fftconv(i1,i2, 'same'))))
% i1=rand(100); i2=rand(100); max(max(abs(1-convn(i1, i2)./fftconv(i1,i2))))
% i1=rand(100,1); i2=rand(100,1); max(abs(1-conv(i1, i2, 'same')./fftconv(i1,i2, 'same')))
