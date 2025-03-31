%% simpleHuygensFresnel1D.m  MN 2024-04-19
% Simple 1D simulation of the far field pattern of an array of point sources
% 
% Usage: [Ez, th, E0, x, ef] = simpleHuygensFresnel1D(x, E0[, option, [value]])
%   Returns:
%     Ez: Complex scalar field amplitude vector at z
%     th: Angle vector corresponding to E values
%
%   Parameters:
%     x: Input axis vector
%     E0: Input scalar electric field vector
%
%     Options:
%       'z', %f: propagation distance (default 100)
%       'th', %f: Specify output angle vector, either as grid, bounds, or span (default pi)
%       'N', %i: Generate grid via linspace(min(th),max(th),N) if numel(th)<3 (default 2^12)
%       'lambda', %f: Specify wavelength (default 1.55e-6)
%       'k', %f: Specify wavenumber (default 2*pi/lambda)
%       'reverse': Treat inputs as far-field plane and outputs as nearfield
%       'elementfactor' | 'ef', [%f]: vector of element factor scaling, corresponding to 'th' grid
%
% TODO:
%   x Demonstrate
%   x Nonuniform input/output size
%   x Reversed operation
%   x Removed multiple-z capability
%   x Change output to V/rad
%   x Element factor
%   x Fix normalization: seems correct
%   x Break up overly large simulations

function [Ez, th, E0, x, ef] = simpleHuygensFresnel1D(x, E0, varargin)
%% Defaults and magic numbers
z = 100;
N = NaN;
th = pi;
lambda = 1.55e-6;     k=2*pi/lambda;
reverse = false;
ef = NaN;


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
        case 'z'
            z = double(nextarg("z"));
            z = z(1);
        case 'th'
            th = double(nextarg("Theta"));
        case 'n'
            N = round(nextarg("N"));
        case 'lambda'
            lambda = double(nextarg("lambda"));
            k=2*pi/lambda;
        case 'k'
            k = double(nextarg("k"));
            lambda=2*pi/k;
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


%% Verify and standardize inputs
if isnan(N); N = 2^12; end

if numel(th) < 2
    th = [-th/2, th/2];
end
if numel(th) == 2
    th = linspace(th(1), th(2), N);
end
th = th(:)';

x = x(:);
E0 = E0(:);
if size(x) ~= size(E0)
    error("x vector and E0 vector size mismatch!");
end

% Flip x and th for reversed simulation
if reverse; [x, th] = deal(th, x); end

% Generate uniform-emission element factor
elementFactor = gradient(th) / mean(gradient(th));

% Apply nonuniform element factor if specified
if ~any(isnan(ef))
    if numel(elementFactor) == numel(ef)
        % Apply with correct shape
        elementFactor = elementFactor .* reshape(ef, size(elementFactor));
    else
        error("Element factor specified with length %i, but mismatch with length %i `th`", numel(elementFactor), numel(ef));
    end
end

% Normalize such that the integral is unity
elementFactor = elementFactor .* trapz(th/2/pi, elementFactor.^2)^-0.5;


%% Basic Huygens-Fresnel principle propagation
% Break down by `th` block as necessary
blockN = numel(th);
doubleBytes = 8;
maxSize = memory().MaxPossibleArrayBytes/doubleBytes;
if maxSize < 8 * numel(E0)*blockN
    blockN = ceil(maxSize / 10 / numel(E0));
end
blocks = 1:numel(th); blocks((end+1):(end+blockN-rem(numel(th), blockN))) = NaN; blocks = reshape(blocks, blockN, []);

% Preallocate
Ez = sparse(numel(th),1);

% Process blocks as needed
for ii = blocks
    ii = ii(~isnan(ii));
    %   Far field is the surface of a circle at z distance from origin
    % rho = sqrt((z^2 * cos(th).^2 + (z * sin(th) - x).^2));    % Straightforward form
    rho = sqrt(x.^2 + z.^2 - 2*x.*z.*sin(th(ii)));  % Numerically simpler form
    hfKernel = exp(1i*k*rho) ./ rho; % Huygens-Fresnel propagation kernel
    Ez(ii) = sqrt(1i/lambda) * sum(E0 .* hfKernel .* elementFactor(ii), 1);
end
Ez = gather(Ez);

% Restore output variables for reversed simulation
if reverse; [x, th] = deal(th, x); end

% Unit transformations
Ez = Ez * z / 2;    % Change to V/rad instead of V/m, including integration over 
                    %   the other angular dimension (4π steradian -> 2π rad)

% P0 = C0*sum(abs(E0).^2);
% Pz = C0*trapz(th, abs(Ez).^2);
% fprintf("Pz/P0 = %.4g / %.4g = %.4f\n", Pz, P0, Pz/P0);


end
