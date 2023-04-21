%% Compute xcorr via FFT
% Michael Nickerson 2019
% Replicates behavior of 'xcorr[n]' via FFT for faster performance
function r = fftxcorr(x, y, maxlag)
    if exist('y', 'var') && numel(y) == 1   % 2nd parameter is 'maxlag'
        maxlag = y;
        clear y;
    end
    if ~exist('y', 'var')
        y = x;
    end
    
    if exist('maxlag', 'var')
        maxlag = (maxlag-1)*(size(x)>1) + 1;  % Make correct dimensions
    else
        maxlag = max(size(x), size(y)) - 1;
    end
    
    n = max([maxlag+(size(x)>1) ; size(x) ; size(y)]);
    n = (2*n-2).*(n>1) + 1;
    
    r = fftshift(ifftn( fftn(x, n) .* conj(fftn(y, n)) ));
    
    % TODO: Impove this logic
    inds = cell(size(n));
    for i = 1:numel(n)
        inds{i} = ceil( max(1, n(i)/2-maxlag(i)) ):ceil( min(n(i), n(i)/2+maxlag(i)) );
    end

    r = r(inds{:});
    
end
%% Tests
% i1=rand(100,1); i2=rand(100,1); max(abs((xcorr(i1,i2)-fftxcorr(i1,i2))./xcorr(i1,i2)))
% i1=rand(100,1); i2=rand(100,1); max(abs((xcorr(i1, i2, 100)-fftxcorr(i1,i2, 100))./xcorr(i1,i2, 100)))
% i1=rand(100,1); i2=rand(100,1); max(abs((xcorr(i1, i2, 50)-fftxcorr(i1,i2, 50))./xcorr(i1,i2, 50)))
% i1=rand(20); i2=rand(20); max(max(abs((xcorr2(i1, i2)-fftxcorr(i1,i2))./xcorr2(i1,i2))))
