%% simpleFraunhofer1D.m  MN 2024-04-26
% Simple 1D simulation of the far field pattern of a linear source
% 
% Usage: [Ez, th] = simpleFraunhofer1D(x, E0[, option, [value]])
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
%
% TODO:
%   x Demonstrate
%   - Compare to HF method
%   - Implement 'reverse' and verify

function [Ez, th, E0, x] = simpleFraunhofer1D(x, E0, varargin)
%% Defaults and magic numbers
N = 2^12;
th = NaN;
lambda = 1.55e-6;
reverse = false;


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
    
    % Look for valid arguments
    switch arg
        case 'th'
            th = double(nextarg("Theta")); th = th(:)';
        case 'n'
            N = round(nextarg("N"));
        case 'lambda'
            lambda = double(nextarg("lambda"));
        case 'k'
            lambda = 2*pi/double(nextarg("k"));
        case {"reverse", "inverse"}
            reverse = true;
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
    E0 = interp1(x, E0, xN, "linear");
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
% Upsample as needed to fully cover desired th
thMax = asin(lambda / dx / 2);

if any(isnan(th)) || (max(th) == thMax && min(th) == thMax && numel(th) == N)
    upN = 1;
else
    upN = max([1 ceil(thMax / max(abs(th)))]);
    if upN > 1; upN = 2*upN; end    % Oversampled for later interpolation
end

% Assemble far field coordinates; v = sin(th)/lambda
thF = asin(lambda * linspace(-0.5,0.5,N*upN) / dx)';

% Edge-pad the near field for FFT
E0 = [zeros(ceil((N*upN - numel(E0))/2), 1); E0; zeros(floor((N*upN - numel(E0))/2), 1)];

% Propagate: just an FFT
if ~reverse
    Ez = ifftshift(fft(fftshift(E0), N*upN));
else
    Ez = fftshift(ifft(ifftshift(E0), N*upN))*N*upN;
end

% Resample to desired size
if upN == 1
    th = thF;
else
    Ez = interp1(thF, Ez, th, "linear");
end


end
