%% simpleHyugensFresnel1D.m  MN 2024-04-19
% Simple 1D simulation of the far field pattern of a linear source
% 
% Usage: [Ez, th, z] = simpleHyugensFresnel1D(x, z, E0[, option, [value]])
%   Returns:
%     Ez: Complex scalar field amplitude vector at z
%     th: Angle vector corresponding to E values
%
%   Parameters:
%     x: Input axis vector; is recentered during simulation
%     z: Propagation distance from center of input; can specify multiple distances
%     E0: Input scalar electric field vector
%
%     Options:
%       'th', double: Specify output angle vector, either as grid, bounds, or span (default pi)
%       'N', %i: Regenerate grid via linspace(min(th),max(th),N) and resample (default 2^12+1)
%       'lambda', %f: Specify wavelength (default 1.55e-6)
%       'k', %f: Specify wavenumber (default 2*pi/lambda)
%
% TODO:
%   x Demonstrate
%   x Nonuniform input/output size

function [Ez, th, z] = simpleHyugensFresnel1D(x, z, E0, varargin)
%% Defaults and magic numbers
N = 2^12+1;
th = pi;
lambda = 1.55e-6;     k=2*pi/lambda;


%% Argument parsing
% Check required inputs
if isempty(x) || ~isa(x, 'double')
    error('Required input "x" is not a double!');
end
if isempty(z) || ~isa(z, 'double')
    error('Required input "z" is not a double!');
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
        case 'th'
            th = double(nextarg("Theta")); th = th(:)';
        case 'n'
            N = round(nextarg("N"));
        case 'lambda'
            lambda = double(nextarg("lambda"));
            k=2*pi/lambda;
        case 'k'
            k = double(nextarg("k"));
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

x = x - mean(x); x = x(:);
E0 = E0(:);
if size(x) ~= size(E0)
    error("x vector and E0 vector size mismatch!");
end

z = z(:)';


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
Ez = NaN(numel(th), numel(z));

for i = 1:numel(z)
    zi = (z(i)^2 * cos(th).^2 + (z(i) * sin(th) - x).^2).^0.5;
    Ez(:, i) = sum(E0 .* exp(1i*k*zi), 1)';
end


end
