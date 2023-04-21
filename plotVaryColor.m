%% Helper function to simplify plotting
% Michael Nickerson 2020
% plotVaryColor(x,y,c[, num, title, xla, yla, zla])
% Standard plot with a varying color along the length
% x, y, c all vectors
function h = plotVaryColor(x, y, c, varargin)
    % Make inputs column vectors
    x = reshape(x,[],1);
    y = reshape(y,[],1);
    c = reshape(c,[],1);
    
    % Generate plot
    if numel(varargin) >= 1; figure(varargin{1}); else; figure(); end
    figureSize(gcf, 800, 650);
    h = surf([x,x], [y,y], zeros(size([x,x])), [c,c], 'EdgeColor', 'interp', 'FaceColor', 'no');
    view(2); shading interp;
    if numel(varargin) >= 5; ylabel(colorbar, varargin{5}, 'FontSize', 16, 'FontName', 'Consolas', 'FontWeight', 'Bold'); zlabel(varargin{5}, 'FontSize', 16, 'FontName', 'Consolas', 'FontWeight', 'Bold'); end
    if numel(varargin) >= 4; ylabel(varargin{4}, 'FontSize', 16, 'FontName', 'Consolas', 'FontWeight', 'Bold'); end
    if numel(varargin) >= 3; xlabel(varargin{3}, 'FontSize', 16, 'FontName', 'Consolas', 'FontWeight', 'Bold'); end
    if numel(varargin) >= 2; title(varargin{2}, 'FontSize', 20, 'FontName', 'Source Sans Pro', 'FontWeight', 'Bold'); end

    drawnow;    % Immediate display
end
