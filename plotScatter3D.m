%% Helper function to simplify plotting
% Michael Nickerson 2019
% plotScatter3D(x,y,z,[c, num, tit, xla, yla, zla, cla])
% Plots a scatterplot
function plotScatter3D(x, varargin)
    % Make x a column vector
    if size(x,2) > size(x,1); x = x'; end;
    
    % If x multidimensional, use that for all possible
    if size(x,2) > 1
        y = x(:,2);
    else
        y = varargin{1}; varargin(1) = [];
    end
    if size(x,2) > 2
        z = x(:,3);
    else
        z = varargin{1}; varargin(1) = [];
    end
    if size(x,2) > 3
        c = x(:,4);
    elseif numel(varargin) >= 1 && numel(varargin{1}) == size(x,1)
        c = varargin{1}; varargin(1) = [];
    else
        c = ones([size(x,1), 1]);
    end
    x = x(:,1);
    
    % First merge duplicate rows
    [~, ii, iDup] = unique([x y z], 'rows');  % Find duplicated rows
    x = x(ii);  y = y(ii);  z = z(ii);        % Drop duplicated rows
    c = accumarray(iDup, c, [], @mean);       % Merge c-values
    
    % Normalize c
    normC = c*sign(mean(c, "omitnan"));
    if numel(unique(normC)) > 2
        normC = normC-min(normC);
    end
    normC = normC./max(normC);
    
    % Generate plot
    if numel(varargin) >= 1; figure(varargin{1}); else; figure(); end;
    figureSize(gcf, 800, 650);
    scatter3(x, y, z, 100, normC, 'filled'); view(3);
    shading interp;
    if numel(varargin) >= 6; ylabel(colorbar, varargin{6}, 'FontSize', 16, 'FontName', 'Consolas', 'FontWeight', 'Bold'); end
    if numel(varargin) >= 5; zlabel(varargin{5}, 'FontSize', 16, 'FontName', 'Consolas', 'FontWeight', 'Bold'); end
    if numel(varargin) >= 4; ylabel(varargin{4}, 'FontSize', 16, 'FontName', 'Consolas', 'FontWeight', 'Bold'); end
    if numel(varargin) >= 3; xlabel(varargin{3}, 'FontSize', 16, 'FontName', 'Consolas', 'FontWeight', 'Bold'); end
    if numel(varargin) >= 2; title(varargin{2}, 'FontSize', 20, 'FontName', 'Source Sans Pro', 'FontWeight', 'Bold'); end
    
    drawnow;    % Immediate display
end
