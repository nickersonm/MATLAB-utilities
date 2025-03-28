%% spgd.m  MN 2024-08-02 -- 2025-02-27
% Implementation of the SPGD algorithm including orthonormal dithering
%  (c) Michael Nickerson 2025
% 
% Usage: [coeff, Jp] = spgd(fhandle, coeff, [, option, [value]])
%   Returns:
%     coeff: final coefficients
%     Jp: final function value, fhandle(coeff)
%
%   Parameters:
%     fhandle: handle to function to evaluate, must accept (coeff) and return scalar double
%     coeff: initial parameter set vector, assumed double
%
%     Options:
%       'clow', [%f]: coefficient lower limit: coeff = max(coeff, clow) (default -inf)
%       'chigh', [%f]: coefficient upper limit: coeff = mod(coeff, chigh) (default inf)
%       'gamma', %f: step size multiple / learning rate (default 0.1)
%       'dither', [%f]: dither weights of coefficients (default 1)
%       'iter', %i: Number of iterations to execute (default inf)
%       'stop', %f: Stopping threshold, fractional change in fhandle(coeff) (default 1e-4)
%       'callback', <handle>: call with (coeff, J) every iteration
%       'quiet': don't show messages
%       'gleak', %f: leakage coefficient (default 1-gamma^2)
%       'min' | 'max': Minimize or maximize (default maximize)
%
% TODO:
%   x Initial implementation
%   x Test
%   x Callback test
%   - Other options?
%       - Cap gradient?

function [coeff, Jp] = spgd(fhandle, coeff, varargin)
%% Defaults and magic numbers

quiet = false;
iter = inf;
stop = 1e-4;
gamma = 0.1;
callback = NaN;
dither = NaN;
clow = -inf;
chigh = inf;
gleak = NaN;
optsign = 1;


%% Argument parsing
% Check required inputs
if isempty(fhandle) || ~isa(fhandle, "function_handle")
    error('Required input "fhandle" is not a function_handle!');
end
if isempty(coeff) || ~isa(coeff, "double")
    error('Required input "coeff" is not a double!');
end
if ~isa(fhandle(coeff), "double") || numel(fhandle(coeff)) > 1
    error('"fhandle" does not return a scalar double!');
end

% Accept a struct.option = value structure
if numel(varargin) > 0 && isstruct(varargin{1})
    coefftruct = varargin{1}; varargin(1) = [];
    varargin = [reshape([fieldnames(coefftruct) struct2cell(coefftruct)]', 1, []), varargin];
end

% Parameter parsing
while ~isempty(varargin)
    arg = lower(varargin{1}); varargin(1) = [];
    if isempty(arg); continue; end
    
    % Look for valid arguments
    switch arg
        case {'iter', 'iterations'}
            iter = int32(abs(round(nextarg('iter'))));
        case {'stop', 'thresh', 'threshold'}
            stop = double(nextarg('stop'));
        case {'gamma', 'gain'}
            gamma = double(nextarg('beta'));
        case {'clow', 'low', 'lowlimit'}
            clow = double(nextarg('coeff low limits'));
        case {'clim', 'high', 'chigh', 'limit', 'limits', 'coefflim'}
            chigh = double(nextarg('coeff high limits'));
        case {'dither', 'weight', 'weights', 'step'}
            dither = double(nextarg('coeff dither weights'));
        case {'gleak', 'leak'}
            gleak = double(nextarg('leakage coefficient'));
        case 'callback'
            callback = nextarg('beta');
            if ~isa(callback, "function_handle")
                error('"callback" is not a valid function_handle');
            end
        case 'max'
            optsign = 1;
        case 'min'
            optsign = -1;
        case 'quiet'
            quiet = true;
        otherwise
            if ~isempty(arg)
                warning('Unexpected option "%s", ignoring', num2str(arg));
            end
    end
end


%% Helper functions, if any
    % Overwrite previous output if passed; quiet if quiet set
    function out = utilDisp(out, varargin)
        if ~quiet
            if numel(varargin) > 0; lastout = varargin{1}; else; lastout = ''; end;
            
            fprintf(repmat('\b', 1, numel(sprintf(lastout))));
            fprintf(out);
        end
    end
    
    % Get the next argument or error
    function arg = nextarg(strExpected)
        if isempty(strExpected); strExpected = ''; end
        if ~isempty(varargin)
            arg = varargin{1}; varargin(1) = [];
        else
            error('Expected next argument "%s", but no more arguments present!', strExpected);
        end
    end
    
    % Modulo both positive and negative
    function u = modbounds(u, clow, chigh)
        for ii = 1:numel(u)
            if u(ii) < clow(ii)
                u(ii) = sign(clow(ii)) * mod(abs(u(ii)), abs(clow(ii)));
            end
            if u(ii) > chigh(ii)
                u(ii) = sign(chigh(ii)) * mod(abs(u(ii)), abs(chigh(ii)));
            end
        end
    end


%% Parameter parsing
iter = int32(max([iter, 1]));
coeff = coeff(:)';
Nu = numel(coeff);

if isnan(gleak); gleak = 1-gamma^2; end

if numel(chigh) < Nu
    chigh(numel(chigh):Nu) = mean(chigh);
end
chigh = chigh(1:Nu);

if numel(clow) < Nu
    clow(numel(clow):Nu) = mean(clow);
end
clow = clow(1:Nu);

if any(isnan(dither)); dither = 1; end
if numel(dither) ~= numel(coeff)
    dither = mean(dither) * ones(size(coeff));
end

function J = cbhandle(coeff)
    J = fhandle(coeff);
    callback(coeff, J);
end
if isa(callback, "function_handle"); fJ = @cbhandle; else; fJ = fhandle; end


%% Generate orthonormal codes
% Using Hadamard codes
codes = hadamard(2^ceil(log2(Nu+1)));
codes = codes(:, 2:Nu+1); % Disregard the first unbalanced codes and drop unnecessary codes

% Alternative: use Gold codes, e.g. https://github.com/gsongsong/matlab-goldcode
% %   Known tap values for n in [5, 6, 7, 9, 10, 11]
% n = [5, 6, 7, 9, 10, 11];
% codes = goldcode(n(max([1, find(log2(Nu) > n, 1)])))';
% [~, ii] = sort(abs(sum(codes, 1)));
% codes = codes(:, sort(ii(1:Nu)));   % Keep most balanced codes

Nh = size(codes,1);


%% Main loop
i=1;
J = ones(Nh, 1) * fJ(coeff);
dJ = diff([J; J(1)]);
x = coeff + gamma * cumsum(codes,1) .* dither;
dx = diff([x(end,:);x],1,1);
imax = floor(intmax/Nh)*Nh;

while i <= iter
    % Define circular indicies
    i0 = mod(i-2, Nh)+1;    % Previous iteration
    i1 = mod(i-1, Nh)+1;    % This iteration
    i2 = mod(i-0, Nh)+1;    % Next iteration
    
    % Apply coefficient for this iteration
    J(i1) = fJ(x(i1,:));
    dJ(i1) = J(i1) - J(i0);
    dx(i1,:) = x(i1,:) - x(i0,:);
    
    % Check for stop
    if abs(dJ(i1)/mean(J)) < stop
        coeff = x(i1,:);
        break;
    end
    
    % Determine best estimate of the gradient
    dJdx = dJ'/dx';
    
    % TODO: cap the gradient?
    % dJdx = sign(dJdx) .* min(abs(dJdx), 1e3);
    
    % Set next iteration's coefficient
    du = gamma * dJdx + codes(i2, :) .* dither;
    u = x(i1,:)*gleak + optsign * du;
    x(i2,:) = modbounds(u, clow, chigh);
    
    i = mod(i, imax) + 1;
end

% Final function value
Jp = fJ(coeff);

end
