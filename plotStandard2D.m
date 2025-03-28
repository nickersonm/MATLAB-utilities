%% plotStandard2D.m  MN 2020-06-29
% Generate a plot with the given data in a standard format
% 
% Usage: [figH, axH, plotH] = plotStandard2D([parameters])
%   Returns:
%     figH:  Handle to figure
%     axH: Handle to axes [x1y1, x2y1, x1y2, x2y2]
%     plotH: Handle to plots
%
%   Parameters:
%     Data: applies to most recently specified plot data only
%       [xv, yv]: x and y column vectors to plot
%       Any unrecognized options are passed to the `plot` command
%       'opts', cell: Include specified cell contents in `plot` command
%       'style', string: Use this style
%       'legend', string: Add a legend <string>
%       'y2': Plot y-data on y2 axis
%       'x2': Plot x-data on x2 axis
%       'samecolor': Keep color the same as the previous plot; not valid for first dataset
%     Figure: applies to all inputs
%       ('figure' | 'fig'), h: Plot using existing or new figure handle h
%       'size', [xSize, ySize]: Make figure this size
%       ('pos' | 'position'), [x y dx dy]: Set figure position directly
%       ('invisible' | 'hide'): Don't show figure
%       'add': Add data to existing plot (hold); requires valid ('figure', h)
%       'logx': Make x-axis logscale
%       'logy': Make y-axis logscale
%       'xlim', [xMin, xMax]: Set x limit (±inf for auto limit)
%       'ylim', [yMin, yMax]: Set y limit (±inf for auto limit)
%       'xlabel', string: Specify x label
%       'ylabel', string: Specify y label
%       'axfont', string: Specify axis font, default 'Consolas'
%       ('axsize' | 'axfs'), double: Specify axis font size, default 14
%       ('axweight' | 'axfw'), string: Specify axis font weight, default 'bold'
%       'title', string: Specify plot title
%       'tfont', string: Specify title font; default 'Source Sans Pro'
%       ('tsize' | 'tfs'), double: Specify title font size; default 20
%       ('tweight' | 'tfw'), string: Specify title font weight; default 'normal'
%       'interpreter', string: Specify title interpreter
%       'legendloc', string: Specify legend location; default 'best'
%       'legendor', ('vertical' | 'horizontal'): Specify legend orientation
%       'x2label', string: specify x2 label
%       'y2label', string: Specify y2 label
%       'logx2': Make x2-axis logscale
%       'logy2': Make y2-axis logscale
%       'x2lim', [x2Min, x2Max]: Set x2 limit (±inf for auto limit)
%       'y2lim', [y2Min, y2Max]: Set y2 limit (±inf for auto limit)
%       'save', filename: Save to filename (detect format from extension)
%       'dpi', dpi: DPI to save with; default 96
%       ('nocolor' | 'noaxcol' | 'norecolor'): Don't recolor axes based on plot colors
%
% TODO:
%   x Symmeterize x2 and y2 behavior/coloring
%   x Check save functionality
%   ~ Synchronize axes when x1y1 combined with x1y2 etc.
%   - Implement varying color along length of line?  Colorbar?
%   - Improve elegance of title assignment code
%   x Fix 'add': currently doesn't work well!  Maybe axes results are random? Fixed using axH.Tag
%   - Add error bar plotting via 'fill' or 'patch'

function [figH, axH, plotH] = plotStandard2D(varargin)
%% Defaults and magic numbers
figH = [];
addData = 0;
plotSize = [800, 600];  % Default
plotPos = NaN;
visible = true;

axLims = ones(4,2) .* [-inf, inf];
axScale = {'linear', 'linear', 'linear', 'linear'};
axRecolor = 1;
xLa = ""; yLa = ""; x2La = ""; y2La = "";
tit = "";
interpreter = "tex";

titF  = 'Source Sans Pro';
titFS = 20;
titFW = 'bold';
axF   = 'Consolas';
axFS  = 14;
axFW  = 'bold';

legendLoc = 'best'; legendOr = 'vertical';
savFile = '';
dpi = 96;

allData = {};
allLegend = {};
allOpts = {};
allAxes = [];
allColor = [];


%% Argument parsing
% Allow passing of cells of options
varargin = flatten(varargin);

% Parameter parsing
while ~isempty(varargin)
    arg = lower(varargin{1}); varargin(1) = [];
    if isempty(arg); continue; end
    
    % Check if argument is additional data
    if isnumeric(arg) && (length(arg) > 1)
        allData{end+1} = checkPlotData(arg);
        allLegend{end+1} = [];
        allOpts{end+1} = {'LineWidth', 3};
        allAxes(end+1,:) = [1, 1];
        allColor{end+1} = [];
        continue;
    end
    
    % Otherwise continue to look for other arguments
    switch arg
        case 'opts'
            allOpts{end} = flatten(varargin{1}, allOpts{end});
            varargin(1) = [];
        case 'style'
            allOpts{end} = flatten(varargin{1}, allOpts{end});
            varargin(1) = [];
        case {'figure', 'fig', 'handle'}
            figH = varargin{1}; varargin(1) = [];
        case 'size'
            plotSize = round(varargin{1}); varargin(1) = [];
        case {'position', 'pos'}
            plotPos = round(varargin{1}, 2); varargin(1) = [];
            if all(plotSize == [800, 600]); plotSize = [NaN, NaN]; end   % Override default size
        case {'invisible', 'hide'}
            visible = false;
        case {'add', 'hold'}
            addData = 1;
        case 'color'
            allColor{end} = double(varargin{1}); varargin(1) = [];
        case 'samecolor'
            allColor{end} = 'same';
        case {'legend', 'name'}
            allLegend{end} = varargin{1}; varargin(1) = [];
        case {'legendloc', 'legendlocation'}
            legendLoc = varargin{1}; varargin(1) = [];
        case {'legendor', 'legendorientation', 'legenddir'}
            legendOr = varargin{1}; varargin(1) = [];
        case {'xlim', 'xlimit', 'xra', 'xrange'}
            axLims(1,:) = double(varargin{1}); varargin(1) = [];
        case {'ylim', 'ylimit', 'yra', 'yrange'}
            axLims(2,:) = double(varargin{1}); varargin(1) = [];
        case {'logx', 'xlog'}
            axScale{1} = 'log';
        case {'logy', 'ylog'}
            axScale{2} = 'log';
        case {'xlabel', 'xla'}
            xLa = varargin{1}; varargin(1) = [];
        case {'ylabel', 'yla', 'y1label'}
            yLa = varargin{1}; varargin(1) = [];
        case {'axisfont', 'axfont', 'axisf', 'axf'}
            axF = varargin{1}; varargin(1) = [];
        case {'axisfontsize', 'axisfs', 'axfs', 'axissize', 'axsize'}
            axFS = double(varargin{1}); varargin(1) = [];
        case {'axisfontweight', 'axisfw', 'axfw', 'axisweight', 'axweight'}
            axFW = string(varargin{1}); varargin(1) = [];
        case {'title', 'tit'}
            tit = varargin{1}; varargin(1) = [];
        case {'titlefont', 'titfont', 'tfont', 'tf'}
            titF = varargin{1}; varargin(1) = [];
        case {'titlefontsize', 'titlefs', 'titfs', 'titlesize', 'titsize', 'tsize', 'tfs'}
            titFS = double(varargin{1}); varargin(1) = [];
        case {'titlefontweight', 'titlefw', 'titfw', 'titleweight', 'titweight', 'tweight', 'tfw'}
            titFW = string(varargin{1}); varargin(1) = [];
        case {'interpreter', 'titletype'}
            interpreter = varargin{1}; varargin(1) = [];
        case {'x2label', 'x2la'}
            x2La = varargin{1}; varargin(1) = [];
        case {'y2label', 'y2la'}
            y2La = varargin{1}; varargin(1) = [];
        case 'y2'
            allAxes(end,2) = 2;
        case 'x2'
            allAxes(end,1) = 2;
        case 'y1'
            allAxes(end,2) = 1;
        case 'x1'
            allAxes(end,1) = 1;
        case {'x2lim', 'x2limit', 'x2ra', 'x2range'}
            axLims(3,:) = double(varargin{1}); varargin(1) = [];
        case {'y2lim', 'y2limit', 'y2ra', 'y2range'}
            axLims(4,:) = double(varargin{1}); varargin(1) = [];
        case {'logx2', 'x2log'}
            axScale{3} = 'log';
        case {'logy2', 'y2log'}
            axScale{4} = 'log';
        case {'save', 'out'}
            savFile = varargin{1}; varargin(1) = [];
        case 'dpi'
            dpi = double(varargin{1}); varargin(1) = [];
        case {'nocolor', 'norecolor', 'noaxcol', 'noaxrecolor'}
            axRecolor = 0;
        otherwise
            if ~isempty(arg)
%                 warning('Unexpected option "%s", adding to plot options of most recent line', num2str(arg));
%                 if isempty(allOpts)
%                     allOpts{end} = arg;
%                 else
                    allOpts{end} = flatten(allOpts{end}, arg);
%                 end
            end
    end
end


%% Dependent parameters or defaults
% Default size if not specified and no existing figure
if any(isnan(plotSize)) && ( isempty(figH) || ~ishandle(figH) || (ishandle(figH) && ~isgraphics(figH)) )
    plotSize = [800,600];
end

% Sort limits
axLims = sort(axLims, 2);


%% Helper functions, if any
    % Regularize plot data
    function plotData = checkPlotData(plotData)
        % Check for transposed input data
        if size(plotData, 2) > 2 && size(plotData, 1) == 2
            plotData = transpose(plotData);
        end
        
        % Check for insufficiently specified plot data
        if size(plotData, 2) < 2
            error('Specified plot data does not contain two columns!');
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
    
    % figureSize to reduce external function dependency and reduce flash if figure exists
    function h = figureSize(h, x, y, visible)
        % Default sizes
        if isempty(x) || isnan(x) || ~isnumeric(x); x = 800; end
        if isempty(y) || isnan(y) || ~isnumeric(y); y = 600; end
        
        % Get real handle if exists
        if isnumeric(h) % Try to find existing numbered figure
            hh = findobj('Number', h, 'Type', 'figure');
            if ~isempty(hh) && isvalid(hh); h = hh; end  % Existing valid figure
        end
        if ~isnumeric(h) && ishandle(h) && isvalid(h)    % Existing valid figure
            h.Visible = visible;
        elseif isnumeric(h) && ~isempty(h)   % New numeric figure
            h = figure(h); h.Visible = visible;
        else    % No or invalid figure specified
            h = figure('visible', visible);
        end
        
        % Get current position and set new size
        pos = get(h, 'Position');
        if ~isnan(x) && ~isnan(y)
            set(h, 'Position', [pos(1) pos(2)+pos(4)-y x y]);
        end
    end


%% Set up plot
% Select/create figure and set size
figH = figureSize(figH, plotSize(1), plotSize(2), visible);

% Set position if specified
if all(~isnan(plotPos)) && numel(plotPos) >= 4 && isnumeric(plotPos)
    figH.Position = plotPos;
end

% Set up axes
if addData ~= 1
    % Clear if not adding data
    clf(figH);
    
    % Build axes: [x1y1, x2y1, x1y2, x2y2], tag with order
    axH = arrayfun(@(i) axes(figH, 'Tag', num2str(i)), 1:4);
    set(axH, 'NextPlot', 'add');
    set(axH(2:end), 'Color', 'none');
else
    % Retreive existing axes and sort to expected order
    axH = figH.Children;
    axH = axH(strcmpi(get(axH, 'type'), 'axes'));
    [~,axI] = sort([axH.Tag]);
    axH = axH(axI);
    
    % Add any additional axes required
    if numel(axH) < 4
        axH = [axH(:), cellfun(@(x) axes(figH), cell(4-numel(axH), 1))];
    end
end

% Set axes locations
set(axH([2,4]), 'XAxisLocation', 'top');
set(axH([3,4]), 'YAxisLocation', 'right');

% Change axes reference to linear array properly addressing axH vector
allAxes = sub2ind([2, 2], allAxes(:,1), allAxes(:,2));


%% Format plot
% Set axis visibility and labels appropriate for this dataset
%   Axis vector is [x1y1, x2y1, x1y2, x2y2]
set(axH([2,3]), 'Visible', 'off');
xlabel(axH(1), xLa, 'FontSize', axFS, 'FontName', axF, 'Interpreter', interpreter, 'FontWeight', axFW);
ylabel(axH(1), yLa, 'FontSize', axFS, 'FontName', axF, 'Interpreter', interpreter, 'FontWEight', axFW);
if any( sum(allAxes == [2,4], 2) ) % X2 exists
    xlabel(axH(4), x2La, 'FontSize', axFS, 'FontName', axF, 'Interpreter', interpreter, 'FontWeight', axFW);
else
    axH(4).XAxis.Visible = 'off';
end
if any( sum(allAxes == [3,4], 2) ) % Y2 exists
    ylabel(axH(4), y2La, 'FontSize', axFS, 'FontName', axF, 'Interpreter', interpreter, 'FontWeight', axFW);
else
    axH(4).YAxis.Visible = 'off';
end

% Scale settings
set(axH([1,3]), 'XScale', axScale{1}); % x1
set(axH([1,2]), 'YScale', axScale{2}); % y1
set(axH([2,4]), 'XScale', axScale{3}); % x2
set(axH([3,4]), 'YScale', axScale{4}); % y2

% Title, making space if needed
%   x2 may get displaced by the title
if (axH(4).XAxis.Visible == "on")
    if ~isempty(tit); title(axH(4), tit, 'FontSize', titFS, 'FontName', titF, 'Interpreter', interpreter, 'FontWeight', titFW); end
    set(axH, 'Position', axH(4).Position);
else
    if ~isempty(tit); title(axH(1), tit, 'FontSize', titFS, 'FontName', titF, 'Interpreter', interpreter, 'FontWeight', titFW); end
    set(axH, 'Position', axH(1).Position);
end

% Link axes scales
% linkprop(axH([1 3]), 'XLim');  % Same X1 axis
% linkprop(axH([4 2]), 'XLim');  % Same X2 axis
% linkprop(axH([1 2]), 'YLim');  % Same Y1 axis
% linkprop(axH([4 3]), 'YLim');  % Same Y2 axis
linkaxes(axH([1 3]), 'x');  % Same X1 axis
linkaxes(axH([4 2]), 'x');  % Same X2 axis
linkaxes(axH([1 2]), 'y');  % Same Y1 axis
linkaxes(axH([4 3]), 'y');  % Same Y2 axis

% Axes limits, if specified
set(axH, 'XLimMode', 'auto'); set(axH, 'YLimMode', 'auto');
set(axH( all([1;0;1;0] * ~isinf(axLims(1,:)), 2) ), 'XLim', axLims(1,:)); % x1
set(axH( all([1;1;0;0] * ~isinf(axLims(2,:)), 2) ), 'YLim', axLims(2,:)); % y1
set(axH( all([0;1;0;1] * ~isinf(axLims(3,:)), 2) ), 'XLim', axLims(3,:)); % x2
set(axH( all([0;0;1;1] * ~isinf(axLims(4,:)), 2) ), 'YLim', axLims(4,:)); % y2


%% Plot data
% Initialize plot handles
plotH = gobjects(size(allData));

% Plot each data vector in turn
for i=1:numel(allData)
    axPlot = axH(allAxes(i));
    
    % Assemble plot command
    plotCmd = flatten(axPlot, allData{i}(:,1), allData{i}(:,2), allOpts{i});
    
    % Update colororder
    axPlot.ColorOrderIndex = max([axH.ColorOrderIndex]);
    
    % Set color if specified
    if ~isempty(allColor{i}) && isnumeric(allColor{i})
        % Normalize if needed
        if any(allColor{i} > 1)
            allColor{i} = allColor{i} / max(allColor{i});
        end
        plotCmd = flatten(plotCmd, 'Color', allColor{i});
    elseif strcmp(allColor{i},'same')
        axPlot.ColorOrderIndex = mod(axPlot.ColorOrderIndex - 2, size(axPlot.ColorOrder, 1))+1;
    end
    
    % Plot
    hold(axPlot, 'on');
    plotH(i) = plot(plotCmd{:});
end


%% Plot-dependent formatting
% Legend, if any specified
if any(~cellfun(@isempty, allLegend))
    % On axis 4 so it's on top
    legend(axH(4), plotH(~cellfun(@isempty, allLegend)), allLegend(~cellfun(@isempty, allLegend)), ...
           'Location', legendLoc, 'FontSize', 16, 'Orientation', legendOr, 'Color', 'white');
end

% Set axes color if exists and recoloring enabled
if (axH(4).XAxis.Visible == "on") && (axRecolor == 1)
    axH(1).XColor = plotH(find( sum(allAxes == [1,3], 2), 1)).Color;
    axH(4).XColor = plotH(find( sum(allAxes == [2,4], 2), 1)).Color;
end
if (axH(4).YAxis.Visible == "on") && (axRecolor == 1)
    axH(1).YColor = plotH(find( sum(allAxes == [1,2], 2), 1)).Color;
    axH(4).YColor = plotH(find( sum(allAxes == [3,4], 2), 1)).Color;
end

% Finish up
hold(axH, 'off'); grid(axH(1), 'on');
drawnow;


%% Save if desired
if ~isempty(savFile)
    if ~(isempty(fileparts(savFile)) || fileparts(savFile) == "") && ~exist(fileparts(savFile), 'dir'); mkdir(fileparts(savFile)); end
    [~, ~, ext] = fileparts(savFile);  ext = char(ext);
    print(figH, savFile, ['-d' ext(2:end)], ['-r' num2str(dpi)]);
end


%% Bring to front if not invisible
if visible
    figure(figH);
end

end
