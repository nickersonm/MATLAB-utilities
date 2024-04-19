%% simplePhasedArray.m  MN 2024-04-19
% Simple 1D simulation of the far field pattern of a phased array
% 
% Usage: [Ez, th, x0, E0] = simpleHyugensFresnel1D(x, ph[, option, [value]])
%   Returns:
%     Ez: Complex scalar field amplitude vector at z
%     th: Angle vector corresponding to E values
%     x0: Constructed near field location vector
%     E0: Constructed near field
%
%   Parameters:
%     x: Locations of uniform emitters
%     ph: scalar dph between emitters or vector exact phase of emitters; can be 2D matrix
%
%     Options:
%       'z', double: Propagation distance from center of input (default 100)
%       'dx', double: Width of emitters (default gradient(x))
%       'th', double: Specify output angle vector as bound, span, or grid (default pi)
%       'N', %i: Simulation linear grid size (default 2^12+1)
%       'lambda', %f: Specify wavelength (default 1.55e-6)
%       'k', %f: Specify wavenumber (default 2*pi/lambda)
%       'plot', handle: Plot results with specified handle
%
% TODO:
%   - Demonstrate

function [Ez, th, x0, E0] = simplePhasedArray(x, ph, varargin)
%% Defaults and magic numbers
N = 2^12+1;
th = pi;
lambda = 1.55e-6;     k=2*pi/lambda;
z = 100;
plotH = []; dx = [];


%% Argument parsing
% Check required inputs
if isempty(x) || ~isa(x, 'double')
    error('Required input "x" is not a double!');
end
if isempty(ph) || ~isa(ph, 'double')
    error('Required input "ph" is not a double!');
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
            z = double(nextarg("z")); z = z(1);
        case 'dx'
            dx = double(nextarg("dx"));
        case 'th'
            th = double(nextarg("Theta")); th = th(:)';
        case 'n'
            N = round(nextarg("N"));
        case 'lambda'
            lambda = double(nextarg("lambda"));
            k=2*pi/lambda;
        case 'k'
            k = double(nextarg("k"));
        case 'plot'
            plotH = nextarg("plot handle");
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
th = th(:);

x = x - mean(x); x = x(:);

if isempty(dx)
    dx = gradient(x);
end
dx = dx(:);
if numel(dx) > 1 && numel(dx) ~= numel(x)
    dx = dx(1);
end

if numel(ph) == 1
    ph = ph * (0:(numel(x)-1))';
end
if ~(size(ph,1) == numel(x)) && (size(ph,2) == numel(x))
    ph = ph';
end
if size(ph,1) ~= numel(x)
    error("Mismatch between 'x' and 'ph' input sizes");
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


%% Build input vectors
x0 = linspace(min(x) - dx(1)/2, max(x) + dx(end)/2, N);
Ex = ((x0 >= (x-dx/2)) & (x0 <= (x+dx/2)));

x0 = x0(:);
Ex = normalize(sum(Ex,1)', "range");


%% Apply phases and compute
Ez = NaN(numel(th), size(ph,2)); E0 = NaN(N, size(ph,2));

for i = 1:size(ph,2)
    dph = sum(((x0' >= (x-dx/2)) & (x0' <= (x+dx/2))) .* ph(:,i), 1)';
    E0(:,i) = Ex .* exp(1i*dph);
    
    [Ez(:,i), th] = simpleHyugensFresnel1D(x0, z, E0(:,i), ...
                                           "N", N, "th", th, "lambda", lambda);
end


%% Plot
if ~isempty(plotH)
    figureSize(plotH, 1200, 500); clf;
    mgn = [0.12, 0.08];
    h = subplot_tight(1,2,1, mgn);
    set(gca, "Position", gca().Position - [0.02 0 0 0]);
    yyaxis left;
    plot(x0, abs(E0).^2, "LineWidth", 2); axis padded;
    grid on;
    xlabel("Emitter Position [m]");
    ylabel("Emitter Amplitude [arb]");
    yyaxis right;
    plot(x0, angle(E0)/pi, ":", "LineWidth", 2); axis padded;
    ylabel("Emitter Phase [\pi rad]");
    title(h, "Near Field", 'FontSize', 14);
    
    h = subplot_tight(1,2,2, mgn);
    set(gca, "Position", gca().Position + [0.02 0 0 0]);
    yyaxis left;
    plot(th/pi, abs(Ez).^2, "LineWidth", 2); axis padded;
    grid on;
    xlabel("Far Field Angle [\pi rad]");
    ylabel("Far Field Amplitude [arb]");
    yyaxis right;
    plot(th/pi, angle(Ez)/pi, ":", "LineWidth", 2); axis padded;
    ylabel("Far Field Phase [\pi rad]");
    title(h, "Far Field", 'FontSize', 14);
    
    h = sgtitle(sprintf('1D OPA Simulation; z = %.4g m', z)); h.FontSize = 16; h.FontWeight = 'bold';
    drawnow;
end


end
