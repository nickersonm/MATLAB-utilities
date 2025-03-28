%% handleField1D.m  MN 2025-03-25
% Expands or collapses appropriate dimensions to handle 1D cases with 2D math
% 
% Requirements:
%   - None
% 
% Usage: [E, x, y, D] = handleField1D(x, y, E[, option[, value]])
%   Returns:
%     x, y: expanded dimension vectors
%     E: expanded field matrix
%     D: expanded dimension
%
%   Parameters:
%     x, y: input dimensions, the first singleton dimension will be expanded
%     E: field (or arbitrary vector) to expand to matrix; dimensions must match
%
%     Options:
%       'collapse': collapses first 2-element dimension
%       'dim', %i: collapse this dimension
%
% TODO:
%   x Implement
%   - Integrate with other functions for test

function [E, x, y, collapseD] = handleField1D(x, y, E, varargin)
%% Defaults and magic numbers
collapseD = NaN;


%% Argument parsing
% Check required inputs
if (isempty(x) || ~isa(x, 'double')) && (isempty(y) || ~isa(y, 'double'))
    error('Neither input "x" and "y" are doubles!');
end
if isempty(E) || ~isa(E, 'double')
    error('Required input "E" is not a double!');
end
if ~all(size(E) == [numel(y), numel(x)])
    error('Required input "E" does not match "x" and "y" size!');
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
        case {'collapse', 'reverse'}
            collapseD = find(size(E) == 2);
        case {'collapsed', 'dim'}
            collapseD = round(nextarg("collapse dimension"));
        otherwise
            if ~isempty(arg)
                warning('Unexpected option "%s", ignoring', num2str(arg));
            end
    end
end

% Dependent processing
if isempty(collapseD)
    error("Could not find dimension to collapse; size(E) = [%s]", num2str(size(E)));
end


%% Helper functions
    % Get the next argument or error
    function arg = nextarg(strExpected)
        if isempty(strExpected); strExpected = ''; end
        if ~isempty(varargin)
            arg = varargin{1}; varargin(1) = [];
        else
            error('Expected next argument "%s", but no more arguments present!', strExpected);
        end
    end


if isnan(collapseD)
%% Determine collapse dimension and expand
    if numel(x) == 1 && numel(y) == 1
        return;
    end
    if numel(x) == 1 && size(E, 2) == 1
        x = [-0.5 0.5] + x;
        E = [E E];
        collapseD = 2;
    elseif numel(y) == 1 && size(E, 1) == 1
        y = [-0.5 0.5] + y;
        E = [E; E];
        collapseD = 1;
    end
else
%% Collapse any expanded dimension
    E = mean(E, collapseD);
    if collapseD == 2
        x = mean(x);
    elseif collapseD == 1
        y = mean(y);
    else
        error("Unknown 'collpaseD'!");
    end
end


end
