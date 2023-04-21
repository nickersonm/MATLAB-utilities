%%
beams = arrayCPA(3)*0.01; z = 1500;
% beams = beams+0.01*rand(size(beams))/5;
beams = [zeros(size(beams,1),1) beams]*rotx(30); beams(:,1) = [];
x = linspace(-0.1,0.1,2^11);  x2 = 5*x;
y = linspace(-0.075,0.075,2^11);  y2 = 5*y;
E0 = efieldGaussianBeam(x, y, beams, 'q', 20i);
Eg = efieldGaussianBeam(x2, y2, beams, 'q', z+20i, 'N', 2^8);

Ez = efieldMeanKernel(x, y, z, E0, 'plot', 2, 'N', 2^8, 'xz', x2, 'yz', y2, 'valcheck', 0);
% 
Ez = efieldDirectFresnel(x, y, z, E0, 'plot', 2, 'N', 2^8, 'xz', x2, 'yz', y2, 'valcheck', 0);
% crop = floor(size(Ez)/8);
% Ez = Ez(crop:end-crop, crop:end-crop); Eg = Eg(crop:end-crop, crop:end-crop);

% disp( angle( Eg(floor(end/2),floor(end/2)) ) / angle( Ez(floor(end/2),floor(end/2)) ) );

%%
figureSize(1, 1200,800);
h = subplot(2,2,1);
surf( abs(Eg).^2 ); shading interp; axis tight; view(2); colorbar;
title(h, 'Gaussian Propagation', 'FontSize', 14);

h = subplot(2,2,2);
surf( abs(Ez).^2 ); shading interp; axis tight; view(2); colorbar;
title(h, 'Numerical Propagation', 'FontSize', 14);

h = subplot(2,2,3);
surf(log10( abs(abs(Ez).^2-abs(Eg).^2)./round(gather(abs(Eg).^2),5) )); shading interp; axis tight; view(3); colorbar;
title(h, 'Log10 Relative Error', 'FontSize', 14);

h = subplot(2,2,4);
surf(log10( abs(abs(Ez).^2-abs(Eg).^2)./max(max(abs(Eg).^2)) )); shading interp; axis tight; view(3); colorbar;
title(h, 'Log10 Normalized Error', 'FontSize', 14);
