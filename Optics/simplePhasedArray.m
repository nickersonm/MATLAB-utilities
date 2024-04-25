%% simplePhasedArray.m  MN 2024-04-19
% Simple 1D simulation of the far field pattern of a phased array
% 
% Usage: [Ez, th, E0, x0] = simplePhasedArray(x, ph[, option, [value]])
%   Returns:
%     Ez: Complex scalar field amplitude vector at z
%     th: Angle vector corresponding to E values
%     E0: Constructed near field
%     x0: Constructed near field location vector
%
%   Parameters:
%     x: Locations of uniform emitters
%     ph: scalar dph between emitters or vector exact phase of emitters; can be 2D matrix
%
%     Options:
%       'P', double: Total power of emitters (default 1)
%       'z', double: Propagation distance from center of input (default 100)
%       'th', double: Specify output angle vector as bound, span, or grid (default pi)
%       'N', %i: Far-field linear grid size (default 2^12+1)
%       'lambda', %f: Specify wavelength (default 1.55e-6)
%       'k', %f: Specify wavenumber (default 2*pi/lambda)
%       'plot', handle: Plot results with specified handle
%       'nonorm': Don't normalize the output by the element count
%
% TODO:
%   x Demonstrate
%   x Normalize input power
%   - Change output to per-radian
%   x Nonuniform input/output size
%   x Don't mesh nearfield - use point emitters

function [Ez, th, E0, x] = simplePhasedArray(x, ph, varargin)
%% Defaults and magic numbers
N = 2^12+1;
P = 1;
th = pi;
lambda = 1.55e-6;
z = 100;
plotH = []; dx = [];
nonorm = false;


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
        case 'p'
            P = double(nextarg("P")); P = P(1);
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
        case 'k'
            lambda = 2*pi/double(nextarg("k"));
        case 'plot'
            plotH = nextarg("plot handle");
        case 'nonorm'
            nonorm = true;
        otherwise
            if ~isempty(arg)
                warning('Unexpected option "%s", ignoring', num2str(arg));
            end
    end
end


%% Verify and standardize inputs
if numel(th) < 2
    th = [-th/2, th/2];
end
if numel(th) == 2
    th = linspace(th(1), th(2), N(end));
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
Ex = 0*x + 1;

% Normalize power
Ex = Ex * (P / sum(gradient(x) .* abs(Ex).^2, "all"))^0.5;


%% Apply phases and compute
Ez = NaN(numel(th), size(ph,2)); E0 = NaN(numel(x), size(ph,2));

for i = 1:size(ph,2)
    E0(:,i) = Ex .* exp(1i*ph(:,i));
    
    Ez(:,i) = simpleHyugensFresnel1D(x, z, E0(:,i), "th", th, "lambda", lambda);
end

% Change units to per radian
%   TODO
% disp(sum(gradient(x0) .* abs(Ex).^2, 1));
if ~nonorm; Ez = Ez / numel(x); end
% disp(sum(gradient(th) .* abs(Ez).^2, 1));


%% Plot
if ~isempty(plotH)
    figureSize(plotH, 1200, 500); clf;
    mgn = [0.12, 0.08];
    h = subplot_tight(1,2,1, mgn);
    set(gca, "Position", gca().Position - [0.02 0 0 0]);
    yyaxis left;
    plot(x, abs(E0).^2, "x", "LineWidth", 2); axis padded;
    grid on;
    xlabel("Emitter Position [m]");
    ylabel("Emitter Amplitude [V/m]");
    yyaxis right;
    plot(x, -angle(E0)/pi, "o", "LineWidth", 2); axis padded;
    ylabel("Emitter Phase [\pi rad]");
    title(h, "Near Field", 'FontSize', 14);
    
    h = subplot_tight(1,2,2, mgn);
    set(gca, "Position", gca().Position + [0.02 0 0 0]);
    yyaxis left;
    plot(th/pi, abs(Ez).^2, "LineWidth", 2); axis padded;
    grid on;
    xlabel("Far Field Angle [\pi rad]");
    ylabel("Far Field Amplitude [V/rad]");
    yyaxis right;
    plot(th/pi, angle(Ez)/pi, ":", "LineWidth", 2); axis padded;
    ylabel("Far Field Phase [\pi rad]");
    title(h, "Far Field", 'FontSize', 14);
    
    h = sgtitle(sprintf('1D OPA Simulation; z = %.4g m', z)); h.FontSize = 16; h.FontWeight = 'bold';
    drawnow;
end


end
