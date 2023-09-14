%% efieldGaussianBeam.m  MN 2018-09-14
% Calculates the E-field intensity for arbitrary gaussian beam sources
% 
% Requirements:
%   - None
% 
% Usage: [E, x, y] = efieldGaussianBeam(x, y, sources[, option, value])
%   Returns:
%     E: Complex field INTENSITY matrix
%     x: x-grid
%     y: y-grid
%
%   Parameters:
%     x, y: Vectors or meshgrids of positions to calculate field
%           If <= 2 elements or 'N' specified, generates grid
%     sources: Mx(3|4) array of source parameters: [x, y, q[, phi, E0]]
%           Optional 'phi' column is per-source phase
%           Optional 'E0' column is per-source relative intensity
%
%     Options:
%       'plot', %i: Plot abs(E)^2 and angle(E) in specified figure
%       'q', %f+1i%f: Complex beam parameter for all sources; allows 
%           2-column 'sources' input
%       'N', %i: [re]generate grid via linspace(min(x),max(x),N)
%       'lambda', %f: Specify wavelength (default 1.5e-6)
%       'k', %f: Specify wavenumber (default 2*pi/lambda)
%
% TODO:
%   - 

function [E, x, y] = efieldGaussianBeam(x, y, sources, varargin)
%% Defaults and magic numbers
figN = NaN;
q = NaN;
N = 2^10;
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
            if argval>0
                figN = round(argval);
            end
        case 'q'
            q = double(argval);
        case 'n'
            N = reshape(round(argval), 1, []);
            x = unique([min(x) max(x)]);
            y = unique([min(y) max(y)]);
        case 'lambda'
            lambda = double(argval);
            k=2*pi/lambda;
        case 'k'
            k = double(argval);
    end
end


%% Helper functions
Erq = @(r, q) conj((1-real(q)/q) * exp(-1i*k*(real(q) + r.^2/(2*q)) )); % Experimentally needs the conujgate
r = @(x,y) sqrt(x.^2+y.^2);


%% Verify and standardize inputs
% Make sure source list is good - assume xy symmetric if only one dimension
% provided
if size(sources,2) == 1
    sources = [sources sources];
end
if size(sources,2) == 2 && ~isnan(q)
    sources = [sources q*ones(size(sources,1),1)];
end
if size(sources,2) < 3
    error('Invalid source size: %i', size(sources));
elseif ~isnan(q)
    sources(:,3) = q;
end
if size(sources,2) < 4
    sources = [sources zeros(size(sources,1),1)]; % Zero per-source phase
end
if size(sources,2) < 5
    sources = [sources ones(size(sources,1),1)]; % Equal starting intensity
end

% Standardize grid inputs
if isscalar(x); x = [-x x]; end     % Change scale to range
if isscalar(y); y = [-y y]'; end
if isvector(x); x = reshape(x, 1, []); end  % Force dimensions if vector
if isvector(y); y = reshape(y, [], 1); end
x = sort(x,2);  % Make monotonic
y = sort(y,1);

% Generate grid if needed
if isscalar(N); N = [N N]; end      % Convert to 2 elements
N(isnan(N)) = mean(N, "omitnan");   % Fill any missing dimensions
if numel(x) == 2
    x = linspace(x(1), x(end), N(1));
end
if numel(y) == 2
    y = linspace(y(1), y(end), N(2))';
end


%% Calculate E fields
E = arrayfun(@(sx,sy,sq,sphi,sI) sI*Erq(r(x-sx, y-sy), sq).*exp(1i*sphi), sources(:,1), sources(:,2), sources(:,3), sources(:,4), sources(:,5), 'UniformOutput', 0);
E = sum(cat(3, E{:}), 3);


%% Optionally plot
if ~isnan(figN)
    figureSize(figN, 1200, 500);
    h = subplot(1,2,1);
    imagesc([min(x(:)) max(x(:))],[min(y(:)) max(y(:))], abs(E).^2); axis image xy; colorbar;
    title(h, 'Intensity', 'FontSize', 14);
    
    h = subplot(1,2,2);
    imagesc([min(x(:)) max(x(:))],[min(y(:)) max(y(:))], angle(E), 'AlphaData', abs(E), 'AlphaDataMapping', 'scaled'); axis image xy;
    colorbar; 
    title(h, 'Phase', 'FontSize', 14);
    
    h = sgtitle('Gaussian Field'); h.FontSize = 14; h.FontWeight = 'bold';
    drawnow;
end


%% Return
x = gather(x);
y = gather(y);
E = gather(E);

end
