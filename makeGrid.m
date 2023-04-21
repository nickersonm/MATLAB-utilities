%% Helper function to simplify plotting
% Michael Nickerson 2019
% [xg, yg, zg] = makeGrid(x,y,z)
% Changed scattered data to gridded data
% x, y, z all vectors
function [x, y, z] = makeGrid(x, y, z)
    % Force column vectors
    if size(x,2) > size(x,1); x = x'; end;
    if size(y,2) > size(y,1); y = y'; end;
    if size(z,2) > size(z,1); z = z'; end;
    
    % First merge duplicate rows
    [~, ii, iDup] = unique([x y], 'rows');  % Find duplicated rows
    x = x(ii);  y = y(ii);                  % Drop duplicated rows
    zi = accumarray(iDup, z, [], @mean);    % Merge z-values
    
    % Get coordinate mapping and reduce x and y
    [x, ~, izx] = unique(x);  [y, ~, izy] = unique(y);
    
    % Make a gridded z and assign correct values
    [~,~,z] = ndgrid(x,y,NaN);  iz = sub2ind(size(z), izx, izy);
    z(iz) = zi;
end
