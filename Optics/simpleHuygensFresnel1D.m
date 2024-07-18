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
%       'N', %i: Generate grid via linspace(min(th),max(th),N) if numel(th)<3 (default 2^12+1)
%       'lambda', %f: Specify wavelength (default 1.55e-6)
%       'k', %f: Specify wavenumber (default 2*pi/lambda)
%       'reverse': Treat inputs as far-field plane and outputs as nearfield
%       'elementfactor' | 'ef', [%f]: vector of element factor scaling, corresponding to 'th' grid
%
% TODO:
%   x Demonstrate
%   x Nonuniform input/output size
%   x Reversed operation
%   - Fix normalization
%       x Removed multiple-z capability
%       x Change output to V/rad
%   - Element factor

function [Ez, th, E0, x, ef] = simpleHuygensFresnel1D(x, E0, varargin)
%% Defaults and magic numbers
z = 100;
N = 2^12+1;
th = pi;
lambda = 1.55e-6;     k=2*pi/lambda;
reverse = false;
ef = [];
C0 = (2*376.73)^-1; % Siemens; C0 == eps0*c/2


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
        case 'z'
            z = double(nextarg("z"));
            z = z(1);
        case 'th'
            th = double(nextarg("Theta")); th = th(:)';
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
            ef = double(nextarg("Element factor")); ef = ef(:)';
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

% TBD: Apply nonuniform element factor if specified
if ~isempty(ef)
    if numel(elementFactor) == numel(ef)
        % Normalize such that the integral is unity
        figure(2); plot(th, ef); hold on;
        ef = ef .* trapz(th/2/pi, ef.^2)^-0.5;
        plot(th, ef); hold off; grid on;
        legend(["Original", "Normalized"]);
        
        % Apply
        elementFactor = elementFactor .* ef;
    else
        error("Element factor specified with length %i, but mismatch with length %i `th`", numel(elementFactor), numel(ef));
    end
end


%% Basic Huygens-Fresnel principle propagation
%   Far field is the surface of a circle at z distance from origin
% rho = sqrt((z^2 * cos(th).^2 + (z * sin(th) - x).^2));    % Straightforward form
rho = sqrt(x.^2 + z.^2 - 2*x.*z.*sin(th));  % Numerically simpler form
hfKernel = exp(-1i*k*rho) ./ rho; % Huygens-Fresnel propagation kernel
Ez = sqrt(1i/lambda) * sum(E0 .* hfKernel .* elementFactor, 1);

% Restore output variables for reversed simulation
if reverse; [x, th] = deal(th, x); end

% Unit transformations
Ez = Ez * z / 2;    % Change to V/rad instead of V/m, including integration over 
                    %   the other angular dimension (4π steradian -> 2π rad)

P0 = C0*sum(abs(E0).^2);
Pz = C0*trapz(th, abs(Ez).^2);
fprintf("Pz/P0 = %.4g / %.4g = %.4f\n", Pz, P0, Pz/P0);


end
