%% Gaussian Propagation test

%% Simulation parameters
z = 500;
% zflat = 20;   % Flatten phase after 10cm
lambda = 1e-6;     k=2*pi/lambda;

NbeamsD = 3;
spacing = 5e-3;
w0 = 2e-3;
steer = 0.0;  % half-waves to steer up or down

N = 2^12;

scale = 1*sqrt(z)*NbeamsD*spacing+6*w0;

x = double(linspace(-scale/2,scale/2,N));
y = x';
% x = x + rand(size(x))*scale1/N/10;   % Add some random sampling to the grid
% y = y + rand(size(y))*scale1/N/10;

figH = figureSize(1, 1200, 1200);


%% Gaussian beam equations
zR = pi*w0^2/lambda;

% % Via physical units
% Rz = @(z) z+zR^2/z;
% wz = @(z) w0*sqrt(1+(z/zR)^2);
% gouy = @(z) atan(z/zR);
% 
% Erz = @(r,z) (w0/wz(z))*exp( -r.^2/wz(z)^2 - 1i*(k*z + k*r.^2/(2*Rz(z)) - gouy(z)) );

% Via complex parameter
q0 = 1i*zR;
qz = @(qi, z) qi + z;
Erq = @(r, q) (1-real(q)/q) * exp(-1i*k*(real(q) + r.^2/(2*q)) );
r = @(x,y) sqrt(x.^2+y.^2);


%% Define sources
% Define beam centers
xcen = (1:NbeamsD); xcen = xcen-mean(xcen);
[xcen, ycen] = meshgrid(xcen);
ikeep = xcen.^2+ycen.^2 <= (NbeamsD/2)^2;
beams = [xcen(ikeep) ycen(ikeep)];

% Summed field function, vectorized version
cell3cat = @(x) cat(3, x{:});
Esum = @(x,y,q) sum(cell3cat(arrayfun(@(bx,by) Erq(r(x+spacing*bx, y+spacing*by), q).*exp(-1i*bx*pi*steer), beams(:,1), beams(:,2), 'UniformOutput', 0)), 3);

% % Alternative - 3x slower
% Esum = @(x,y,q) reshape(sum(Erq(r(bsxfun(@plus, x(:), spacing*beams(:,1)'), bsxfun(@plus, y(:), spacing*beams(:,2)')), q).*exp(-1i*beams(:,2)*pi*steer)', 2), size(x));
%
% % Unrolled version for clarity
% function E = Esum(x,y,q)
%     E = zeros(N,N);
%     for ii=1:size(beams,1)
%         E = E + Erq(r(x+spacing*beams(ii,1), y+spacing*beams(ii,2)), q)*exp(-1i*beams(ii,2)*pi*steer);
%     end
% end


%% Plot near field
q = q0;

% Add beam fields
u1 = Esum(x,y,q);

% Plot
subplot(2,2,1);
% surf(x,y,abs(u1).^2); axis image; view(2); shading interp;
imagesc(x([1 end]),y([1 end]), abs(u1).^2); axis image xy;
colorbar; drawnow;

subplot(2,2,3);
% surf(x,y,angle(u1), 'AlphaData', abs(u1), 'FaceAlpha', 'interp'); axis image; view(2); shading interp;
imagesc(x([1 end]),y([1 end]), angle(u1), 'AlphaData', abs(u1), 'AlphaDataMapping', 'scaled'); axis image xy;
colorbar; drawnow;


%% Flatten phase after zflat
% q0 = q0 + (1/(1i*imag(1/(q0+zflat))) - (q0+zflat));


%% Far field
q = q0 + z;

% Add beam fields
u2 = Esum(x,y,q);

figure(figH);
subplot(2,2,2);
imagesc(x([1 end]),y([1 end]), abs(u2.^2)); axis image xy;
colorbar;

subplot(2,2,4);
imagesc(x([1 end]),y([1 end]), angle(u2), 'AlphaData', abs(u2), 'AlphaDataMapping', 'scaled'); axis image xy;
colorbar; drawnow;


%% Attempt at slices
% Nslice = 20;
% 
% figure(2); figureSize(gcf, 1200, 800);
% for zz = 0:z/Nslice:z
%     q = q0 + zz;
%     E = Esum(x,y,q);
%     colorVal = angle(E);
% %     colorVal = -abs(E).^2;
% %     alphaVal = (abs(E).^0.5/max(max(abs(E).^0.5)));
%     alphaVal = abs(E).^2/max(max(abs(E).^2));
%     surf(ones(size(x))*zz, x, y, colorVal, 'AlphaData', alphaVal, 'FaceAlpha', 'interp' ); shading interp; axis tight; hold on;
%     drawnow;
% end
% % colorbar;
% hold off;

