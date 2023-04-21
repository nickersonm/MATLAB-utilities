%% Helper function to simplify plotting
% Michael Nickerson 2019
% plotHeatSurf(x,y,z[, num, title, xla, yla, zla])
% Plots a surface with 2D heatmap-type view
% x, y, z all vectors
function h = plotHeatSurf(x, varargin)
    % Make x a column vector
    if size(x,2) > size(x,1); x = x'; end
    
    % If x multidimensional, use that for all possible
    if size(x,2) > 1 && size(x,2) < 4
        y = x(:,2);
    else
        y = varargin{1}; varargin(1) = [];
    end
    if size(x,2) > 2 && size(x,2) < 4
        z = x(:,3);
    else
        z = varargin{1}; varargin(1) = [];
    end
    x = x(:,1);
    
    % Make a z grid if not already in the right format
    if size(z, 2) < 2
        % First merge duplicate rows
        [~, ii, iDup] = unique([x y], 'rows');  % Find duplicated rows
        x = x(ii);  y = y(ii);                  % Drop duplicated rows
        zz = accumarray(iDup, z, [], @mean);    % Merge z-values

        % Get coordinate mapping and reduce x and y
        [x, ~, izx] = unique(x);  [y, ~, izy] = unique(y);

        % Make a gridded z and assign correct values
        [~,~,z] = ndgrid(x,y,NaN);  iz = sub2ind(size(z), izx, izy);
        z(iz) = zz;
    end

    % Generate plot
    if numel(varargin) >= 1; figure(varargin{1}); else; figure(); end
    figureSize(gcf, 800, 650);
    h = surf(x, y, z'); view(2); axis vis3d;
    shading interp; %colormap jet;
    if numel(varargin) >= 5; ylabel(colorbar, varargin{5}, 'FontSize', 16, 'FontName', 'Consolas', 'FontWeight', 'Bold'); zlabel(varargin{5}, 'FontSize', 16, 'FontName', 'Consolas', 'FontWeight', 'Bold'); end
    if numel(varargin) >= 4; ylabel(varargin{4}, 'FontSize', 16, 'FontName', 'Consolas', 'FontWeight', 'Bold'); end
    if numel(varargin) >= 3; xlabel(varargin{3}, 'FontSize', 16, 'FontName', 'Consolas', 'FontWeight', 'Bold'); end
    if numel(varargin) >= 2; title(varargin{2}, 'FontSize', 20, 'FontName', 'Source Sans Pro', 'FontWeight', 'Bold'); end

    drawnow;    % Immediate display
end
