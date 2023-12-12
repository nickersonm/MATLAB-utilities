%% efieldMeanKernel.m  MN 2018-09-14
% Directly propagates a given scalar E-field via Mean-Kernel matricies
%   Implements method from https://doi.org/10.1364/OE.23.026853 (Section 3)
%   This essentially treats each grid-box as a rect source and sums the
%       propagated fields; similar to EfieldGaussianBeam
%   Slightly modified to treat field intensity instead of field, to allow
%       for resolution and scale independence
%   Slightly slower than direct Fresnel FFT, but higher accuracy
% 
% Requirements:
%   - None
% 
% Usage: [Ez, xz, yz, Ei, x, y] = efieldMeanKernel(x, y, z, Ei[, option, value])
%   Returns:
%     Ez: Complex field INTENSITY matrix at z
%     xz: x-grid at z
%     yz: y-grid at z
%     Ei: utilized input field
%     x: utilized input x-grid
%     y: utilized input y-grid
%
%   Parameters:
%     x, y: Input axes: range, vector, or meshgrid
%         If size doesn't match Ei, generates based on range
%     z: Scalar distance to propagate
%     Ei: Initial complex field INTENSITY matrix
%
%     Options:
%       'plot', %i: Plot abs(E)^2 and angle(E) in specified figure
%       'plotphase', 0: Don't plot phase
%       'N', %i: Output grid size; accepts 2 dimensions for nonuniform grid
%       'lambda', %f: Specify wavelength (default 1.5e-6)
%       'k', %f: Specify wavenumber (default 2*pi/lambda)
%       'xout' | 'x2' | 'xz', [%f]: Specify Ez x axis; range, vector, or meshgrid
%       'yout' | 'y2' | 'yz', [%f]: Specify Ez y axis; range, vector, or meshgrid
%       'valcheck', 1: Reverse-propagate to recover original Ei as well
%           (requires plot to be useful)
%
% TODO:
%   x Custom output grid
%   x Nonuniform grid
%   - Speed up
%   - Fix field vs. intensity confusion
%   - Auto-split into reasonably sized subarrays when input exceeds memory
%   - Implement general EMA (section 2) for arbitrary ABCD transforms?

function [Ez, xz, yz, Ei, x, y] = efieldMeanKernel(x, y, z, Ei, varargin)
%% Defaults and magic numbers
figN = NaN;
N = 0;
lambda = 1.5e-6;     k=2*pi/lambda;
xz = []; yz = [];
valcheck = 0;
plotphase = 1;


%% Argument parsing
% Accept a struct.option = value structure
if numel(varargin) > 0 && isstruct(varargin{1})
    paramStruct = varargin{1}; varargin(1) = [];
    varargin = [reshape([fieldnames(paramStruct) struct2cell(paramStruct)]', 1, []), varargin];
end

if mod(numel(varargin),2)   % I always use "'flag', value" even for boolean commands
    error('Odd number of optional inputs!');
end
% Optional alterations
for i = 1:2:length(varargin)
    arg = lower(varargin{i});
    argval = varargin{i+1};
    switch arg
        case {'plot', 'figure', 'fig'}
            if argval > 0
                figN = round(argval);
            end
        case {'plotphase', 'phase'}
            plotphase = argval;
        case 'n'
            N = reshape(round(argval), 1, []);
        case 'lambda'
            lambda = double(argval);
            k=2*pi/lambda;
        case 'k'
            k = double(argval);
        case {'xout', 'x2', 'xz', 'xf'}
            xz = double(argval);
        case {'yout', 'y2', 'yz', 'yf'}
            yz = double(argval);
        case 'valcheck'
            valcheck = argval > 0;
    end
end


%% Helper functions
    function x = stdgrid(x, dim, sz)
        % x = range, vector, or meshgrid
        % dim = dimension for vector
        % sz = full size to compare
        
        % Select direction for reshape operations
        dims = circshift({1,[]}, dim);
        
        % Change scale to range
        if isscalar(x); x = [-x x]; end
        
        % Force dimensions if vector
        if isvector(x); x = reshape(x, dims{:}); end
        
        % Generate grid if not properly sized vector or meshgrid
        if ~((size(x, dim) == sz(dim)) || (numel(x) == prod(sz)))
            x = reshape(linspace(min(x(:)), max(x(:)), sz(dim)), dims{:});
        end
        
        x = double(x);
    end


%% Verify and standardize inputs
% Standardize input grid
x = stdgrid(x, 2, size(Ei));
y = stdgrid(y, 1, size(Ei));

% If Ei is too small, add a buffer and regenerate xy
if any(size(Ei) <= 3)
    E0 = zeros(size(Ei)+4);
    x = x*size(E0,1)/size(Ei,1)/2.0;    % Rescale to new boundaries
    y = y*size(E0,1)/size(Ei,1)/2.0;    % Rescale to new boundaries
    
    E0(3:end-2,3:end-2) = Ei;  Ei = E0;     % Inelegant
    x = stdgrid(x, 2, size(Ei));
    y = stdgrid(y, 1, size(Ei));
end

% Standardize output grid orientation and range
if isvector(xz); xz = reshape(xz, 1, []); end   % Correct orientation
if isvector(yz); yz = reshape(yz, [], 1); end
if isempty(xz); xz = [min(x(:)) max(x(:))]; end % Use input range if not specified
if isempty(yz); yz = [min(y(:)) max(y(:))]; end

% Check output grid size
if isscalar(N); N = [N N]; end      % Convert to 2 elements
N(isnan(N)) = 0;                    % Next steps rely on 0 == generate

N = [size(xz,2) size(yz,1)] .* [size(xz,2)>2 size(yz,1)>2] .* (N==0) + N;  % Use xz,yz size if specified and not a range
N = size(Ei) .* (N==0) + N;                 % Use input size as fallback

% Standardize or generate output grid
xz = stdgrid(xz, 2, N);
yz = stdgrid(yz, 1, N);


%% Run calculations
H0 = sqrt(exp(1i*k*z)/(1i*lambda*z));

% % Direct calculation
% dx = (gradient(x)' .* gradient(xz)).^0.5;
% dy = (gradient(y) .* gradient(yz)').^0.5;
% x2x1 = (xz - x');
% y2y1 = (yz' - y);
% Hx = H0 * dx .* exp( (1i*pi/(lambda*z)) * x2x1.^2 ) .* sinc(x2x1 .* dx / (lambda*z));
% Hy = H0 * dy .* exp( (1i*pi/(lambda*z)) * y2y1.^2 ) .* sinc(y2y1 .* dy / (lambda*z));

% Method to reduce numerical precision errors:
dxsq = ( gradient(x)' .* ones(size(xz)) ).^2;
x2x1sqz = (xz - x').^2 /z;
Hx = H0 * dxsq.^0.5 .* exp( (1i*pi/(lambda)) * x2x1sqz ) .* sinc( (z*x2x1sqz.*dxsq).^0.5 );

dysq = ( gradient(y) .* ones(size(yz')) ).^2;
y2y1sqz = (yz' - y).^2 /z;
Hy = H0 * dysq.^0.5 .* exp( (1i*pi/(lambda)) * y2y1sqz ) .* sinc( (z*y2y1sqz.*dysq).^0.5 );

% Resultant E-field
Ez = transpose(Hy) * Ei * Hx;   % Note: ' is conjugate transpose!

% Validity check
if valcheck == 1
    % TODO: This doesn't scale properly due to field vs intensity mismatch;
    % should eventually figure out solution!
%     Ei = Ei - conj(Hy) * Ez * Hx';
    Ei = Ei - efieldMeanKernel(xz, yz, -z, Ez, 'xz', x, 'yz', y);
end


%% Optionally plot
if ~isnan(figN)
    mgn = [0.10, 0.05];
    m = double(plotphase ~= 0) + 1;
    
    figureSize(figN, 1200, 800/2*m + 100); clf(figN);
    
    h = subplot_tight(m,2,1, mgn);
    surf(x,y, abs(Ei).^2); shading flat; axis tight; view(2); colorbar;
    title(h, 'Ei Intensity', 'FontSize', 14);
    
    h = subplot_tight(m,2,2, mgn);
    surf(xz,yz, abs(Ez).^2); shading flat; axis tight; view(2); colorbar;
    title(h, 'Ez Intensity', 'FontSize', 14);
    
    if m > 1
        h = subplot_tight(m,2,3, mgn);
        surf(x,y, angle(Ei), 'AlphaData', abs(Ei), 'AlphaDataMapping', 'scaled', 'FaceAlpha', 'flat'); shading flat; axis tight; view(2); colorbar;
        title(h, 'Ei Phase', 'FontSize', 14);
        
        h = subplot_tight(m,2,4, mgn);
        surf(xz,yz, angle(Ez), 'AlphaData', abs(Ez), 'AlphaDataMapping', 'scaled', 'FaceAlpha', 'flat'); shading flat; axis tight; view(2); colorbar;
        title(h, 'Ez Phase', 'FontSize', 14);
    end
    
    sgtitle(sprintf('Mean-Kernel Propagation; z = %.4g', z), 'FontSize', 16, 'FontWeight', 'bold');
    drawnow;
end


%% Return
xz = gather(xz);
yz = gather(yz);
Ez = gather(Ez);

end
