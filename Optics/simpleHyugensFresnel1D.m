%% simpleHyugensFresnel1D.m  MN 2024-04-19
% Simple 1D simulation of the far field pattern of a linear source
% 
% Usage: [Ez, th, z] = simpleHyugensFresnel1D(x, E0[, option, [value]])
%   Returns:
%     Ez: Complex scalar field amplitude vector at z
%     th: Angle vector corresponding to E values
%
%   Parameters:
%     x: Input axis vector; is recentered during simulation
%     E0: Input scalar electric field vector
%
%     Options:
%       'z', double: propagation distance (default 100)
%       'th', double: Specify output angle vector, either as grid, bounds, or span (default pi)
%       'N', %i: Regenerate grid via linspace(min(th),max(th),N) and resample (default 2^12+1)
%       'lambda', %f: Specify wavelength (default 1.55e-6)
%       'k', %f: Specify wavenumber (default 2*pi/lambda)
%       'reverse': Treat inputs as far-field plane and outputs as nearfield
%
% TODO:
%   x Demonstrate
%   x Nonuniform input/output size
%   x Reversed operation

function [Ez, th, E0, x] = simpleHyugensFresnel1D(x, E0, varargin)
%% Defaults and magic numbers
z = 100;
N = 2^12+1;
th = pi;
lambda = 1.55e-6;     k=2*pi/lambda;
reverse = false;


%% Argument parsing
% Check required inputs
if isempty(x) || ~isa(x, 'double')
    error('Required input "x" is not a double!');
end
if isempty(E0) || ~isa(E0, 'double')
    error('Required input "E0" is not a double!');
end

% Allow passing of cells of options
varargin = flatten(varargin);

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
        case 'th'
            th = double(nextarg("Theta")); th = th(:)';
        case 'n'
            N = round(nextarg("N"));
        case 'lambda'
            lambda = double(nextarg("lambda"));
            k=2*pi/lambda;
        case 'k'
            k = double(nextarg("k"));
        case "reverse"
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
th = th(:)';

x = x(:);
E0 = E0(:);
if size(x) ~= size(E0)
    error("x vector and E0 vector size mismatch!");
end

z = z(:)';

% Flip x and th for reversed simulation
if reverse; [x, th] = deal(th, x); end


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
    
    % Flatten a nested cell
    function flatCell = flatten(varargin)
        flatCell = {};
        for j=1:numel(varargin)
            if iscell(varargin{j})
                flatCell = [flatCell flatten(varargin{j}{:})];
            else
                flatCell = [flatCell varargin(j)];
            end
        end
        flatCell = flatCell( ~cellfun(@isempty, flatCell) );
    end


%% Basic Hyugens-Fresnel principle
Ez = NaN(numel(th)*(1-reverse) + numel(x)*reverse, numel(z));

for i = 1:numel(z)
    zi = (z(i)^2 * cos(th).^2 + (z(i) * sin(th) - x).^2).^0.5;
    Ez(:, i) = sum(E0 .* exp(1i*k*zi), 1)';
end

% Restore output variables for reversed simulation
if reverse; [x, th] = deal(th, x); end


end
