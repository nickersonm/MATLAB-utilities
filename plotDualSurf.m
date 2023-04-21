%% Helper function to simplify plotting
% Michael Nickerson 2019
% plotDualSurf(x,y,z, num, title, xla, yla, zla)
% Plots a surface with 2D and 3D view
function plotDualSurf(x,y,z, num, title, xla, yla, zla)
    % Generate Plot with two views
    figure(num); figureSize(gcf, 1200, 500);
    subplot(1,2,1); surf(x, y, z); view(2); axis vis3d;
    shading interp; colormap jet; ylabel(colorbar, zla);
    xlabel(xla); ylabel(yla); zlabel(zla);

    subplot(1,2,2); surf(x, y, z); view(3); axis vis3d;
    shading interp; colormap jet;
    xlabel(xla); ylabel(yla); zlabel(zla);

    % Build title axes and title.
    axes('Position', [0, 1, 1, 0.05], 'Visible', 'off' ) ;
    text(0.5, -0.25, title, 'FontSize', 14', 'FontWeight', 'Bold', 'VerticalAlignment', 'Top', 'HorizontalAlignment', 'Center' ) ;

    drawnow;    % Immediate display
end