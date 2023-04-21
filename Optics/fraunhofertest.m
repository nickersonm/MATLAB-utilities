%% Fraunhofer diffraction test

%% Set up simulation
z = 1;
lambda = 1;

NbeamsD = 4;
spacing = 1;
w0 = 0.3;

N = 2^11;
Nfft = min([N*8 2^14]);

scale = NbeamsD*spacing+6*w0;
[x,y] = meshgrid(linspace(-scale/2,scale/2,N));
res = scale/N;

figure(1); figureSize(gcf, 1200, 550);


%% Build near field
gauss = @(x,y) exp(-(x.^2+y.^2)/w0^2);

% Define beam centers
xcen = (1:NbeamsD); xcen = xcen-mean(xcen);
[xcen, ycen] = meshgrid(xcen);
ikeep = xcen.^2+ycen.^2 <= (NbeamsD/2)^2;
beams = spacing*[xcen(ikeep) ycen(ikeep)];

% Add beam fields - could not figure out easy way to vectorize
u1 = zeros(N,N);
for ii=1:size(beams,1)
    u1 = u1 + gauss(x+beams(ii,1), y+beams(ii,2))*exp(-1i*beams(ii,2)*pi);
end

% Plot
subplot(1,2,1);
surf(x,y,abs(u1).^2); axis equal; view(2); shading interp; drawnow;


%% Far field
% FFT at increased resolution, then remove irrelevant parts
iEnd = [1:N/2 Nfft-N/2+1:Nfft];
u2 = fftn(u1,[Nfft Nfft])/Nfft; u2 = u2(iEnd, iEnd);
u2 = fftshift(u2).*(exp(1i*z*2*pi/lambda)/(1i*lambda*z)).*exp(1i*pi.*(x.^2+y.^2)/(lambda*z));

% Fourier coordinates
[x,y] = meshgrid( linspace(-0.5,0.5,N)*N/(Nfft*res) );

% Scale constants
x = x*lambda*z; y = y*lambda*z;

subplot(1,2,2);
surf(x, y,abs(u2.^2)); axis equal; view(2); shading interp;
xlim([-2 2]); ylim([-2 2]);
