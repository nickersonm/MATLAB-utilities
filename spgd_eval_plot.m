%% spgd_eval_plot.m  MN 2025-02-27
% Plot SPGD evaluation
%  (c) Michael Nickerson 2025
% 
% Usage: spgd_eval_plot(coeff, J, [, option, [value]])
%   Parameters:
%     coeff: coefficient value(s)
%     J: function value
%
%     Options:
%       'figure', %i: plot to specified figure handle (default 1)
%
% TODO:
%   x Initial implementation
%   x Test

function spgd_eval_plot(coeff, J, varargin)
%% Defaults and magic numbers
figN = 1;


%% Argument parsing
% Check required inputs
if isempty(coeff) || ~isa(coeff, "double")
    error('Required input "coeff" is not a double!');
end
if isempty(J) || ~isa(J, "double")
    error('Required input "J" is not a double!');
end

% Accept a struct.option = value structure
if numel(varargin) > 0 && isstruct(varargin{1})
    coefftruct = varargin{1}; varargin(1) = [];
    varargin = [reshape([fieldnames(coefftruct) struct2cell(coefftruct)]', 1, []), varargin];
end

% Parameter parsing
while ~isempty(varargin)
    arg = lower(varargin{1}); varargin(1) = [];
    
    % Look for valid arguments
    switch arg
        case {'fig', 'figure'}
            figN = nextarg('figure number');
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


%% 2-axis plot
h = figure(figN);

% Plot coeffs
yyaxis left;
grid on;
ax = gca();

if numel(ax.Children) ~= numel(coeff) || any(cell2mat(cellfun(@size, {ax.Children.XData}, 'UniformOutput', false)) ~= cell2mat(cellfun(@size, {ax.Children.YData}, 'UniformOutput', false)))
    clf; yyaxis left;
    hold on;
    for i = 1:numel(coeff)
        plot(1, coeff(i));
    end
    hold off;
else
    for i=1:numel(coeff)
        [ax.Children(i).XData(end+1), ax.Children(i).YData(end+1)] = ...
            deal(ax.Children(i).XData(end) + 1, coeff(i));
    end
end

% Plot J
yyaxis right;
ax = gca();
if isempty(ax.Children)
    plot(1, J);
else
    [ax.Children(1).XData(end+1), ax.Children(1).YData(end+1)] = ...
        deal(ax.Children(1).XData(end) + 1, J);
end

drawnow;

end
