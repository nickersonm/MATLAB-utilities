%% Helper function to simplify plotting
% Michael Nickerson 2019
% ax = figureTitle(h, titleText, topMargin, titOpts{:})
% Adds a new invisible axis with custom text to figure handle h
%   TODO: Allow axes to be passed instead
function ax = figureTitle(h, titleText, topMargin, varargin)
%     if ~isempty(varargin); varargin = varargin{1}; end
    if isempty(h)
        h = gcf;
    end
    if ~exist('topMargin', 'var') || isempty(topMargin)
        topMargin = 0.05;
    end
    
    % Get real handle if needed
    if isnumeric(h)
        hh = findobj('Number', h, 'Type', 'figure');
        
        % If not valid, make a new figure
        if ~isempty(hh) && isvalid(hh)
            h = hh;
        end
    end
    h = figure(h);    
    
    % Make new axes; TODO: get existing axis?
    ax = axes(h, 'Position', [0, 1, 1, topMargin], 'Visible', 'off') ;
    
    % Add text
    text(ax, 0.5, -0.5, titleText, 'FontSize', 20, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', 'FontName', 'Source Sans Pro', varargin{:} ) ;
    
    drawnow;    % Immediate display
end