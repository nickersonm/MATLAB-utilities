%% Fraunhofer diffraction test

%% Simulation parameters
% Physical
z = 10;
lambda = 1e-6;     k=2*pi/lambda;

% Source
NbeamsD = 4;
spacing = 0.5e-3;
w0 = .2e-3;
steer = 1.0;  % half-waves to steer up or down

% Grid
N = 2^11;
scale1 = 2*z*(NbeamsD*spacing+6*w0);
scale2 = NbeamsD*spacing+6*w0;

Escale = (scale1/N)^2;


figH = figureSize(3, 1200, 800);


%% Build near field
% Define beam sources
q0 = gaussianQ(lambda, w0);
beams = arrayCPA(NbeamsD);
beamPhi = steer*pi*beams(:,2);
beams = [beams.*spacing q0*ones(size(beams,1),1)];

% Get initial E field
[x1, y1] = meshgrid(linspace(-scale1/2,scale1/2,N));
% x1 = x1 + rand(size(x1))*scale1/N/10;   % Add some random sampling to the grid
% y1 = y1 + rand(size(y1))*scale1/N/10;

E0 = efieldGaussianBeam(x1, y1, [beams beamPhi], 'k', k, 'q', q0);


% Plot initial field
figure(figH); subplot(2,2,1);
imagesc(x1([1 end]),y1([1 end]), abs(E0).^2); axis image xy;
colorbar; drawnow;

subplot(2,2,3);
imagesc(x1([1 end]),y1([1 end]), angle(E0), 'AlphaData', abs(E0), 'AlphaDataMapping', 'scaled'); axis image xy;
colorbar; drawnow;


%% Transformations
E1 = E0;
% % Attempt to flatten phase at 2m
% zflat = 0; zapply = 0.1;
% Eflat = efieldGaussianBeam(x1, y1, beams, 'k', k, 'q', q0+zflat);
% E10 = efieldGaussianBeam(x1, y1, beams, 'k', k, 'q', q0+zapply);
% 
% r = (x1.^2+y1.^2+(zapply-zflat)^2).^0.5;
% hfKernel = -1i * exp(1i*k*r) / (lambda*(zapply-zflat));
% phAp = angle(fftconv(abs(Eflat), conj(hfKernel), 'same')) .* ((x1.^2 + y1.^2) < (scale1/2/z)^2);
% 
% % Plot transformation field
% figureSize(2, 1200, 800);
% subplot(2,2,3);
% imagesc(x1([1 end]),y1([1 end]), phAp); axis image xy; colorbar;
% 
% subplot(2,2,4);
% imagesc(x1([1 end]),y1([1 end]), abs(Eflat).^2); axis image xy; colorbar; drawnow;
% 
% E1 = efieldGaussianBeam(x1, y1, [beams beamPhi], 'k', k, 'q', q0+zapply) .* exp(-1i*phAp);
% z = z-zapply;
% 
% subplot(2,2,1);
% imagesc(x1([1 end]),y1([1 end]), abs(E1).^2); axis image xy; colorbar;
% 
% subplot(2,2,2);
% imagesc(x1([1 end]),y1([1 end]), angle(E1)-angle(E10), 'AlphaData', abs(E1), 'AlphaDataMapping', 'scaled'); axis image xy;
% colorbar; drawnow;


%% Far field via convolution
% % Convolution kernels
% % Huygens or Fresnel-Kirchhoff PSF
% r = (x1.^2+y1.^2+z^2).^0.5;
% % hfKernel = -1i * exp(1i*k*r) ./ (lambda*r); % Full Huygens
% hfKernel = -1i * exp(1i*k*r) / (lambda*z);  % Partial paraxial simplification in denominator
% E2 = fftconv(E1, conj(hfKernel), 'same') * Escale;  % Conjugate required to match gaussian propagation
% 
% % % Alternative kernels / PSF
% % hKernel = -1i * exp(1i*k*z) / (lambda*z) * exp(1i*k*(x1.^2+y1.^2)/(2*z)); % Fresnel 
% % E2h = fftconv(E1, hKernel, 'same') * Escale;
% % figure; imagesc((abs(E2-E2h))./abs(E2)); colorbar;

E2 = efieldDirectFresnel(unique(x1), unique(y1), z, E1, 'k', k);

figure(figH); subplot(2,2,2);
imagesc(x1([1 end])/z,y1([1 end])/z, abs(E2.^2)); axis image xy; colorbar;

subplot(2,2,4);
imagesc(x1([1 end])/z,y1([1 end])/z, angle(E2), 'AlphaData', abs(E2), 'AlphaDataMapping', 'scaled'); axis image xy; colorbar; drawnow;


%% Compare to exact gaussian solution
Egauss = efieldGaussianBeam(x1, y1, [beams beamPhi], 'k', k, 'q', q0+z);
figure(figureSize(2, 1000, 800));
subplot(2,2,1); imagesc(x1([1 end]),y1([1 end]), abs(Egauss.^2) ); axis image xy; colorbar;
subplot(2,2,2); imagesc(x1([1 end]),y1([1 end]), (abs(E2).^2-abs(Egauss).^2)/max(max(abs(E2).^2)) ); axis image xy; colorbar;

subplot(2,2,3); imagesc(x1([1 end]),y1([1 end]), angle(Egauss), 'AlphaData', abs(Egauss), 'AlphaDataMapping', 'scaled'); axis image xy; colorbar; drawnow;
subplot(2,2,4); imagesc(x1([1 end]),y1([1 end]), angle(E2)-angle(Egauss) ); axis image xy; colorbar;
