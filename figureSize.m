%% Helper function to simplify plotting
% Michael Nickerson 2019
% figureSize(h, x, y)
% Changes the size of figure h without changing top-left position
function h = figureSize(h, varargin)
    x=[]; y=[]; visible=[];
    % Process inputs
    while numel(varargin) > 0
        if isempty(x) && isnumeric(varargin{1})
            x = varargin{1};
            if numel(x) > 1; y = x(2); x = x(1); end
        elseif isempty(y) && isnumeric(varargin{1})
            y = varargin{1};
        elseif isempty(visible)  && (isnumeric(varargin{1}) || islogical(varargin{1}))
            visible = varargin{1};
        elseif strcmp(varargin{1}, 'hide') || strcmp(varargin{1}, 'invisible')
            visible = false;
        end
        varargin(1) = [];
    end
    
    % Defaults
    if isempty(x) || isnan(x) || ~isnumeric(x); x = 800; end
    if isempty(y) || isnan(y) || ~isnumeric(y); y = x*0.75; end
    if isempty(visible); visible = true; end
    
    % Get real handle if exists
    if isnumeric(h) % Try to find existing numbered figure
        hh = findobj('Number', h, 'Type', 'figure');
        if ~isempty(hh) && isvalid(hh); h = hh; end  % Existing valid figure
    end
    if ~isnumeric(h) && ishandle(h) && isvalid(h)    % Existing valid figure
        h = figure(h); h.Visible = visible;
    elseif ~isempty(h) && isnumeric(h)   % New numeric figure
        h = figure(h); h.Visible = visible;
    else    % No or invalid figure specified
        h = figure('visible', visible);
    end
    
    % Get current position and set new size
    pos = get(h, 'Position');
    set(h, 'Position', [pos(1) pos(2)+pos(4)-y x y]);
end