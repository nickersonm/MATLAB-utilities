%% filename.m
%   (c) Michael Nickerson YYYY-MM-DD
% Brief description of function
% 
% Requirements:
%   - List of required files
%   - And required state of equipment
% 
% Usage: [outputs, if, any] = functionName(required, inputs[, option, [value]])
%   Returns:
%     output: Description
%
%   Parameters:
%     required: Description of requirement
%
%     Options:
%       'flagname'[, valueOrType]: Description of functionality (default value)
%
% TODO:
%   - Any known TODOs

function [dataOut1, dataOut2] = functionTemplateMN(requiredIn, varargin)
%% Defaults and magic numbers

quiet = 0;
myParam = [];


%% Argument parsing
% Check required inputs
if isempty(requiredIn) || ~isa(requiredIn, 'double')
    error('Required input "requiredIn" is not a double!');
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
        case {"singleopt", "alias"}
            myParam = double(nextarg("my required parameter"));
        case "quiet"
            quiet = 1;
        otherwise
            if ~isempty(arg)
                warning("Unexpected option '%s', ignoring", num2str(arg));
            end
    end
end


%% Dependent variable processing, if any


%% Helper functions, if any
    % Overwrite previous output if passed; quiet if quiet set
    function out = utilDisp(out, varargin)
        if ~quiet
            if numel(varargin) > 0; lastout = varargin{1}; else; lastout = ''; end;
            
            fprintf(repmat('\b', 1, strlength(lastout)));
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


%% Process stuff
dataOut1 = mean(requiredIn);


%% Analyze stuff
dataOut2 = 1./dataout1;


end
