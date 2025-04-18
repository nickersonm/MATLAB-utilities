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
%       'P', %f: Total power of emitters (default 1)
%       'E0', [%f]: Electric field of individual emitters, must match numel(x)
%       'z', %f: Propagation distance from center of input (default 100)
%       'th', %f: Specify output angle vector as bound, span, or grid (default pi)
%       'N', %i: Far-field linear grid size (default 2^12+1)
%       'lambda', %f: Specify wavelength (default 1.55e-6)
%       'k', %f: Specify wavenumber (default 2*pi/lambda)
%       'plot', handle: Plot results with specified handle
%       'nocenter': Don't recenter x vector
%       'elementfactor' | 'ef', [%f]: vector of element factor scaling, corresponding to 'th' grid
%       'norm': Normalize output by comparing to entire 2pi emission
%       'fft' | 'fraunhofer': Use Fraunhofer propagation (FFT) instead of Huygens-Fresnel; ignores z
%
% TODO:
%   x Demonstrate
%   x Normalize input power
%   x Nonuniform input/output size
%   x Don't mesh nearfield - use point emitters
%   x Don't recenter input offsets
%   x Change output to V/rad (via simpleHuygensFresnel1D)
%   x Element factor (via simpleHuygensFresnel1D)
%   x Fix normalization: simply normalize to total 2pi emission
%       - TBD: real physical normalization?
%   x Add simpleFaunhofer1D option ("fft") => normalization is different

function [Ez, th, E0, x] = simplePhasedArray(x, ph, varargin)
%% Defaults and magic numbers
N = 2^12+1;
P = 1;
th = pi;
lambda = 1.55e-6;
z = 100;
plotH = NaN;
nocenter = false;
C0 = (2*376.73)^-1; % Siemens; C0 == eps0*c/2
ef = NaN;
norm = false;
fraunhofer = false;
E0 = NaN;


%% Argument parsing
% Check required inputs
if isempty(x) || ~isa(x, 'double')
    error('Required input "x" is not a double!');
end
if isempty(ph) || ~isa(ph, 'double')
    error('Required input "ph" is not a double!');
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
        case "p"
            P = double(nextarg("P")); P = P(1);
        case "z"
            z = double(nextarg("z")); z = z(1);
        case {"e0", "e", "ei"}
            E0 = double(nextarg("E0"));
            if numel(E0) > 1 && numel(E0) ~= numel(x)
                error("Number of 'E0' elements does not match number of 'x' elements!");
            end
        case "th"
            th = double(nextarg("Theta")); th = th(:)';
        case "n"
            N = round(nextarg("N"));
        case "lambda"
            lambda = double(nextarg("lambda"));
        case "k"
            lambda = 2*pi/double(nextarg("k"));
        case "plot"
            plotH = nextarg("plot handle");
        case "nocenter"
            nocenter = true;
        case {"elementfactor", "element", "ef"}
            ef = double(nextarg("Element factor")); ef = ef(:)';
        case {"norm", "normalize"}
            norm = true;
        case {"fft", "fraunhofer"}
            fraunhofer = true;
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
    th = linspace(th(1), th(2), N(end));
end
th = th(:);

x = x(:);
if ~nocenter; x = x - mean(x); end

if numel(ph) == 1
    ph = ph * (0:(numel(x)-1))';
end
if ~(size(ph,1) == numel(x)) && (size(ph,2) == numel(x))
    ph = ph';
end
if size(ph,1) ~= numel(x)
    error("Mismatch between 'x' and 'ph' input sizes");
end

% Propagation function
if fraunhofer
    prop1D = @(E0, th, ef) simpleFraunhofer1D(x, E0, "th", th, "ef", ef, "lambda", lambda);
else
    prop1D = @(E0, th, ef) simpleHuygensFresnel1D(x, E0, "th", th, "ef", ef, "z", z, "lambda", lambda);
end


%% Build input vectors
if isnan(E0); E0 = 1; end
Ex = 0*x + E0(:);

% Normalize power
%   Total power = C0 * sum(abs(Ex).^2)
Ex = Ex .* (P / (C0 * sum(abs(Ex).^2)))^0.5;


%% Apply phases and compute
Ez = NaN(numel(th), size(ph,2)); E0 = NaN(numel(x), size(ph,2));

for i = 1:size(ph,2)
    E0(:,i) = Ex .* exp(1i*ph(:,i));
    
    Ez(:,i) = prop1D(E0(:,i), th, ef);
    
    if norm
            warning("Normalizing");
            % Normalize by comparing to entire 2pi emission
            th2 = linspace(-pi, pi, ceil(2*pi/mean(diff(th))));
            if ~any(isnan(ef))
                ef2 = interp1(th, ef, th2, "linear", 0);
            else
                ef2 = ef;
            end
            [Ez0, th0] = prop1D(E0(:,i), th2, ef2);
            Pz0 = C0 * trapz(th0, abs(Ez0).^2);
            Ez(:,i) = Ez(:,i) * (P / Pz0)^0.5;
    end
end

% P0 = C0 * sum(abs(Ex).^2);
% Pz = C0 * trapz(th, abs(Ez).^2);
% fprintf("Pz/P0 = %.4g / %.4g = %.4f\n", Pz, P0, Pz/P0);


%% Plot
if ~isempty(plotH) && ~any(isnan(plotH))
    figureSize(plotH, 1200, 500); clf;
    mgn = [0.12, 0.08];
    h = subplot_tight(1,2,1, mgn);
    set(gca, "Position", gca().Position - [0.02 0 0 0]);
    yyaxis left;
    plot(x, abs(E0), "x", "LineWidth", 2); axis padded;
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
    plot(th/pi, abs(Ez), "LineWidth", 2); axis padded;
    grid on;
    xlabel("Far Field Angle [\pi rad]");
    ylabel("Far Field Amplitude [V/rad]");
    yyaxis right;
    plot(th/pi, angle(Ez)/pi, ":", "LineWidth", 1); axis padded;
    ylabel("Far Field Phase [\pi rad]");
    title(h, "Far Field", 'FontSize', 14);
    
    h = sgtitle(sprintf('1D OPA Simulation; z = %.4g m', z)); h.FontSize = 16; h.FontWeight = 'bold';
    drawnow;
end


end
