%% simpleFraunhofer1D.m  MN 2024-04-26
% Simple 1D simulation of the far field pattern of a linear source
% 
% Usage: [Ez, th, E0, x] = simpleFraunhofer1D(x, E0[, option, [value]])
%   Returns:
%     Ez: Complex scalar field amplitude vector at z
%     th: Angle vector corresponding to E values
%
%   Parameters:
%     x: Input axis vector; if not uniform, will resample with count N
%     E0: Input scalar electric field vector
%
%     Options:
%       'th', double: Specify output angle vector, either as grid, bounds, or span (default determined by x)
%       'N', %i: Regenerate grid via linspace(min(th),max(th),N) and resample (default 2^12)
%       'lambda', %f: Specify wavelength (default 1.55e-6)
%       'k', %f: Specify wavenumber (default 2*pi/lambda)
%       'reverse': Treat inputs as far-field plane and outputs as nearfield
%       'elementfactor' | 'ef', [%f]: vector of element factor scaling, corresponding to 'th' grid
%
% TODO:
%   x Demonstrate
%   x Compare to HF method
%   x Implement 'reverse' and verify
%   x Element factor
%   - Normalization: not quite correct, can't figure out
%   x Nearfield subsampling to increase farfield extent
%       x Implement and test
%       - Verify 'reverse' still works with this

function [Ez, th, E0, x] = simpleFraunhofer1D(x, E0, varargin)
%% Defaults and magic numbers
N = 2^12;
th = NaN;
lambda = 1.55e-6;
reverse = false;
ef = NaN;
resample = false;


%% Argument parsing
% Check required inputs
if isempty(x) || ~isa(x, 'double')
    error('Required input "x" is not a double!');
end
if isempty(E0) || ~isa(E0, 'double')
    error('Required input "E0" is not a double!');
end

% Accept a struct.option = value structure
if numel(varargin) > 0 && isstruct(varargin{1})
    paramStruct = varargin{1}; varargin(1) = [];
    varargin = [reshape([fieldnames(paramStruct) struct2cell(paramStruct)]', 1, []), varargin];
end

% Parameter parsing
while ~isempty(varargin)
    arg = lower(varargin{1}); varargin(1) = [];
    if isempty(arg); continue; end
    
    % Look for valid arguments
    switch arg
        case 'th'
            th = double(nextarg("Theta"));
        case 'n'
            N = round(nextarg("N"));
        case 'lambda'
            lambda = double(nextarg("lambda"));
        case 'k'
            lambda = 2*pi/double(nextarg("k"));
        case {"reverse", "inverse"}
            reverse = true;
        case {"elementfactor", "element", "ef"}
            ef = double(nextarg("Element factor"));
            if isnan(N); N = numel(ef); end
        otherwise
            if ~isempty(arg)
                warning('Unexpected option "%s", ignoring', num2str(arg));
            end
    end
end

% Verify and standardize inputs
if numel(th) < 2
    th = [-th/2, th/2];
end
if max(abs(th)) > pi/2
    error("Fraunhofer propagation cannot address greater than ±π/2");
end
if numel(th) == 2
    th = linspace(th(1), th(2), N);
end
N = numel(th);
th = th(:);

x = x(:);
E0 = E0(:);
if size(x) ~= size(E0)
    error("x vector and E0 vector size mismatch!");
end

% Resample input if not a uniform grid
dx = mean(diff(x));
if (std(diff(x))/dx > 1e-3) && (std(diff(sin(x)))/sin(dx) > 1e-3)
    xN = linspace(min(x), max(x), max([N, numel(x)]))';
    E0 = interp1(x, E0, xN, "nearest");
    x = xN(:); clear("xN");
end


%% Helper functions, if any
    % Get the next argument or error
    function arg = nextarg(strExpected)
        if isempty(strExpected); strExpected = ''; end
        if ~isempty(varargin)
            arg = varargin{1}; varargin(1) = [];
        else
            error('Expected next argument "%s", but no more arguments present!', strExpected);
        end
    end


%% Fraunhofer propagation: FFT
% Check coordinate bounds
thMax = asin(lambda / dx / 2);

% Subsample nearfield as needed to fully cover desired farfield
if thMax < max(abs(th))
    % Use integer multiples of original coordinates only for ideal subsampling
    upN = ceil(max(abs(th)) / thMax);
    xN = linspace(min(x), max(x), numel(x)*upN)';
    E0 = interp1(x, E0, xN, "nearest");
    x = xN(:); clear("xN");
    
    % Redetermine thMax
    dx = dx / upN;
    thMax = real(asin(lambda / dx / 2));
end

% Determine far-field coordinates
if any(isnan(th)) || (max(th) == thMax && min(th) == thMax && numel(th) == N)
    % Not specified or bounds match exactly
    upN = 1;
else
    % Other cases: make sure the minimum resolution is met
    upN = max([1 ceil(thMax/N / min(diff(th)))]);
    if upN > 1 && upN < 5; upN = 2*upN; end    % Oversampled for later interpolation
    resample = true;    % Only looking at a subset of the desired angle here; resample for return
end

% Assemble far field coordinates; v = sin(th)/lambda
thF = asin(lambda * linspace(-0.5,0.5,N*upN) / dx)';

% Remove evanescent propagation angles from calculations
iProp = thF==real(thF);
thF(~iProp) = NaN;

% Edge-pad the near field for FFT
E0 = [zeros(ceil((N*upN - numel(E0))/2), 1); E0; zeros(floor((N*upN - numel(E0))/2), 1)];

% Generate uniform-emission element factor
elementFactor = gradient(thF) / mean(gradient(thF(iProp)));
elementFactor([find(iProp,1,'first'), find(iProp,1,'last')]) = 0;

% Apply nonuniform element factor if specified
if ~any(isnan(ef))
    if numel(elementFactor) == numel(ef)
        % Apply with correct shape
        elementFactor = elementFactor .* interp1(th, ef, thF, "linear", 0);
    else
        error("Element factor specified with length %i, but mismatch with length %i `th`", numel(elementFactor), numel(ef));
    end
end

% Normalize such that the integral is unity
elementFactor(iProp) = elementFactor(iProp) .* trapz(thF(iProp)/2/pi, elementFactor(iProp).^2)^-0.5;

% Propagate: just an FFT
if ~reverse
    Ez = elementFactor .* ifftshift(fft(fftshift(E0), N*upN)) / dx^0.5;
else
    Ez = fftshift(ifft(ifftshift(elementFactor .* E0), N*upN))*N*upN / dx^0.5;
end
% figure(3); plotyy(thF, abs(Ez), thF, angle(Ez)); axis tight; drawnow;

% Reduce to non-evanescent portion
thF = thF(iProp);
Ez = Ez(iProp);

% Resample to desired size
if upN == 1 && ~resample
    th = thF;
else
    Ez = interp1(thF, Ez, th, "linear", 0);
end


end
