%%
lambda = 1.55e-6; k = 2*pi/lambda;
z = 1000; z0 = 750e-6;   % 1mm is practicable for powerphotonics
MFD = 2e-6;
spacing = 5*MFD;
beams = arrayCPA(4);
% beams = beams.^2.*sign(beams)*0.5/max(beams(:))^2+beams;
% beams = beams+rand(size(beams))/5;
% beams = [zeros(size(beams,1),1) beams]*rotx(30); beams(:,1) = [];
% beams = [beams 1*ones(size(beams,1),2)*pi.*beams(:,2)];
beams(:,1:2) = spacing*beams(:,1:2);

x0 = 5*spacing*((size(beams,1)^0.5)/2+2); y0 = x0;

[E0, x0, y0] = efieldGaussianBeam(x0, y0, beams, 'q', gaussianQ(lambda, MFD, z0), 'N', 2^9, 'plot', 0);
totalP = sum((gradient(y0)*gradient(x0)) .* abs(E0).^2, 'all');
ph0 = -angle(E0); %ph0 = smoothn(unwrap_phase(gather(ph0)), 10);    % Perfect microlens at z=z0
E0 = E0 .* exp(1i*ph0);

% E0 = E0 .* exp(-1i*(k*(x0.^2+y0.^2)/(1.0*2.0*z0)));   % Spherical lens
% TODO: ideal lens is phi = (2*pi/lambda) * (f + phimax*lambda/(2*pi) - sqrt(x.^2 + f^2))

% xz = linspace(-1,1,2^8); xz = xz.^2 .* sign(xz);
xz = 50;
yz = xz;

[Ez, xz, yz] = efieldMeanKernel(x0, y0, z-z0, E0, 'plot', 2, 'xz', xz, 'yz', yz, 'N', 2^9);

% Ez0 = efieldGaussianBeam(xz, yz, beams, 'q', gaussianQ(1.55e-6, MFD, z), 'N', 2^9, 'plot', 3);

centerR = 0.5/2;
centerI = xz.^2 + yz.^2 < centerR^2;
centerP = mean(abs(Ez(centerI)).^2, 'all')*pi*centerR^2;
% fprintf('Center %.4g urad power fraction: %.4g\n', centerR*2*1000, centerP/totalP);

return;

%% Reduce sidelobes with phase mask / microlens
zflat = 10e-3;
zsc = 2;

Eflat = efieldGaussianBeam(zsc*x0, zsc*y0, beams./zsc, 'q', gaussianQ(lambda, MFD, zflat), 'N', 2^9);
Eflat = Eflat .* exp(-1i*angle(Eflat));
E0 = efieldMeanKernel(zsc*x0, zsc*y0, z0-zflat, real(Eflat), 'plot', 2, 'xz', x0, 'yz', y0);
ph0 = -angle(E0);
figure(1); surf( ph0*lambda/(2*pi*1.55) ); shading interp; axis tight; colorbar;

E0 = efieldGaussianBeam(x0, y0, beams, 'q', gaussianQ(lambda, MFD, z0), 'N', 2^9);
E0 = E0 .* exp(-1i*ph0);
efieldMeanKernel(x0, y0, zflat-z0, E0, 'plot', 3, 'xz', zsc*x0, 'yz', zsc*y0, 'N', 2^9);


%% Spherical fit ph0
ph0 = smoothn(unwrap_phase(gather(ph0)), 10);
sagprofile = 1e6*ph0(256,:)/k; sagprofile = sagprofile - max(sagprofile) - 2;

ft = fittype( 'a - b*(x^2+1)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'Bisquare';
opts.StartPoint = [1 50];

% Fit model to data.
sagfit = fit( 1e3*gather(x0'), sagprofile', ft, opts );
sagfit = feval(sagfit, 1e3*x0)';

figure(1); subplot(1,2,1);
plot(1e3*x0,sagprofile, 1e3*x0,sagfit,':', 'LineWidth', 2);
xlabel('X (mm)'); ylabel('Sag (um optical pathlength)'); title('Freeform Lens');
subplot(1,2,2);
plot(1e3*x0,sagprofile-sagfit, 'LineWidth', 2);
xlabel('X (mm)'); ylabel('Deviation (um optical pathlength)'); title('Deviation from Spherical');


%% Test phase unwraps
% phiz = angle(Ez); phiz((abs(Ez)/max(abs(Ez),[],'all')) <1e-5) = NaN;
phiz = ph0;

figureSize(1, 1200, 400);
phizU = unwrap_phase(gather(phiz));
h = subplot(1,2,1);
surf(phizU);
shading flat; axis tight; view(2); colorbar;
title(h, 'Unwrapped Phase', 'FontSize', 14); drawnow;

h = subplot(1,2,2);
surf(angle(exp(1i .* phizU))-phiz );
shading flat; axis tight; view(2); colorbar;
title(h, 'Unwrapping Error', 'FontSize', 14); drawnow;


%%
% figureSize(1, 1200,800);
% h = subplot(2,2,1);
% surf( abs(Eg).^2 ); shading interp; axis tight; view(2); colorbar;
% title(h, 'Gaussian Propagation', 'FontSize', 14);
% 
% h = subplot(2,2,2);
% surf( abs(Ez).^2 ); shading interp; axis tight; view(2); colorbar;
% title(h, 'Numerical Propagation', 'FontSize', 14);
% 
% h = subplot(2,2,3);
% surf(log10( abs(abs(Ez).^2-abs(Eg).^2)./round(gather(abs(Eg).^2),5) )); shading interp; axis tight; view(3); colorbar;
% title(h, 'Log10 Relative Error', 'FontSize', 14);
% 
% h = subplot(2,2,4);
% surf(log10( abs(abs(Ez).^2-abs(Eg).^2)./max(max(abs(Eg).^2)) )); shading interp; axis tight; view(3); colorbar;
% title(h, 'Log10 Normalized Error', 'FontSize', 14);
