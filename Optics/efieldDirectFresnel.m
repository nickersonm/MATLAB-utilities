%% efieldDirectFresnel.m  MN 2018-09-14
% Directly propagates a given scalar E-field via Fresnel convolution
% 
% Requirements:
%   - None
% 
% Usage: [Ez, xz, yz] = efieldDirectFresnel(x, y, z, Ei[, option, value])
%   Returns:
%     Ez: Complex field amplitude matrix at z
%     xz: x-grid at z
%     yz: y-grid at z
%
%   Parameters:
%     x, y: Vectors or meshgrids of positions to calculate field
%           If <= 2 elements, generates grid with same size as Ei
%     z: Scalar distance to propagate
%     Ei: Matrix of initial field
%
%     Options:
%       'plot', %i: Plot abs(E)^2 and angle(E) in specified figure
%       'N', %i: Regenerate grid via linspace(min(x),max(x),N) and resample
%           Also accepts 2 dimensions for nonuniform grid
%       'lambda', %f: Specify wavelength (default 1.5e-6)
%       'k', %f: Specify wavenumber (default 2*pi/lambda)
%
% TODO:
%   - 

function [Ez, x, y] = efieldDirectFresnel(x, y, z, Ei, varargin)
%% Defaults and magic numbers
figN = NaN;
N = NaN;
lambda = 1.5e-6;     k=2*pi/lambda;


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
        case {'plot', 'figure'}
            figN = round(argval);
        case 'n'
            N = round(argval);
        case 'lambda'
            lambda = double(argval);
            k=2*pi/lambda;
        case 'k'
            k = double(argval);
    end
end


%% Verify and standardize inputs

% Generate grid if needed
if isscalar(x); x = [-x x]; end
if isscalar(y); y = [-y y]; end
if numel(x) < size(Ei,2)
    x = linspace(min(x(:)), max(x(:)), size(Ei,2));
end
if numel(y) < size(Ei,1)
    y = linspace(min(y(:)), max(y(:)), size(Ei,1));
end

% Force dimensions if vector
if isvector(x)
    x = reshape(x, 1, []);
end
if isvector(y)
    y = reshape(y, [], 1);
end

% Make monotonic
[x, xi] = sort(x,2);
[y, yi] = sort(y,1);
Ei = Ei(yi, xi);

% If new N specified, regenerate grid and resample Ei
if isscalar(N); N = [N N]; end
if ~isnan(N)
    xz = double(linspace(min(x(:)), max(x(:)), N(1)) );
    yz = double(linspace(min(y(:)), max(y(:)), N(2)) )';
    
    Ei = interp2(x, y, Ei, xz, yz, 'linear', 0);
    x = xz; y = yz; clear xz yz;
end


%% Run calculations
Escale = (max(x(:))-min(x(:))) * (max(y(:))-min(y(:)))/numel(Ei);

hfKernel = -1i * exp(1i*k*(x.^2+y.^2+z^2).^0.5) / (lambda*z);  % Huygens or Fresnel-Kirchhoff PSF, partial parabolic approximation
x = gather([min(x(:)) max(x(:))]);
y = gather([min(y(:)) max(y(:))]);

Ez = fftconv(Ei, hfKernel, 'same') * Escale;


%% Optionally plot
if ~isnan(figN)
    figureSize(figN, 1200, 800);
    h = subplot(2,2,1);
    imagesc(x,y, abs(Ei).^2); axis image xy; colorbar;
    title(h, 'Ei Intensity', 'FontSize', 14);
    
    h = subplot(2,2,3);
    imagesc(x,y, angle(Ei), 'AlphaData', abs(Ei), 'AlphaDataMapping', 'scaled'); axis image xy; colorbar;
    title(h, 'Ei Phase', 'FontSize', 14);
    
    h = subplot(2,2,2);
    imagesc(x,y, abs(Ez).^2); axis image xy; colorbar;
    title(h, 'Ez Intensity', 'FontSize', 14);
    
    h = subplot(2,2,4);
    imagesc(x,y, angle(Ez), 'AlphaData', abs(Ez), 'AlphaDataMapping', 'scaled'); axis image xy; colorbar;
    title(h, 'Ez Phase', 'FontSize', 14);
    
    h = sgtitle(sprintf('Fresnel Propagation; z = %.4g', z)); h.FontSize = 14; h.FontWeight = 'bold';
    drawnow;
end


%% Return
x = gather(x);
y = gather(y);
Ez = gather(Ez);

end
