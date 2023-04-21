clear;

%%
N = 2^11;
dWG = 0.02;
d = 0.2;
z = 1500;


%%
beams = dWG*arrayCPA(3);
x = d*linspace(-1,1,N); y = x;
E0 = efieldGaussianBeam(x, y, beams, 'q', 20i);

Ezg  = efieldGaussianBeam(x, y, beams, 'q', z+20i);
EzMK = efieldMeanKernel(x, y, z, E0);
EzDF = efieldDirectFresnel(x, y, z, E0);


%%
figureSize(1, 1200,800);
h = subplot(2,2,1);
surf( abs(EzDF).^2 ); shading interp; axis tight; view(2); colorbar;
title(h, 'Fresnel FFT Propagation', 'FontSize', 14);

h = subplot(2,2,2);
surf( abs(EzMK).^2 ); shading interp; axis tight; view(2); colorbar;
title(h, 'MeanKernel Propagation', 'FontSize', 14);

h = subplot(2,2,3);
surf(log10( abs(abs(EzDF).^2-abs(Ezg).^2)./max(max(abs(Ezg).^2)) )); shading interp; axis tight; view(3); colorbar;
title(h, 'Fresnel FFT: Log10 Normalized Error', 'FontSize', 14);

h = subplot(2,2,4);
surf(log10( abs(abs(EzMK).^2-abs(Ezg).^2)./max(max(abs(Ezg).^2)) )); shading interp; axis tight; view(3); colorbar;
title(h, 'MeanKernel: Log10 Normalized Error', 'FontSize', 14);
