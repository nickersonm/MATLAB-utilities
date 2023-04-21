%% Setup
lambda = 1.55e-6; k = 2*pi/lambda;
z = 1000; z0 = 0e-6;   % 1mm is practicable for powerphotonics
xz = 10;  yz = xz;


%% Test simple spacing
spacings = 3:0.1:8;
D = zeros(size(spacings));  Pf = D;

for ii = 1:numel(spacings)
    MFD = 5e-6;
    spacing = spacings(ii)*MFD;
    beams = arrayCPA(4);
    beams(:,1:2) = spacing*beams(:,1:2);
    x0 = 2*spacing*((size(beams,1)^0.5)/2+2); y0 = x0;

    [E0, x0, y0] = efieldGaussianBeam(x0, y0, beams, 'q', gaussianQ(lambda, MFD, z0), 'N', 2^9, 'plot', 0);
    totalP = sum((gradient(y0)*gradient(x0)) .* abs(E0).^2, 'all');

    [Ez, xz, yz] = efieldMeanKernel(x0, y0, z, E0, 'xz', xz, 'yz', yz, 'N', 2^8, 'valcheck', 0);

    centerR = 0.5/2;
    centerI = xz.^2 + yz.^2 < centerR^2;
    centerP = mean(abs(Ez(centerI)).^2, 'all')*pi*centerR^2;
    
    D(ii) = 2*max(beams(:,1));
    Pf(ii) = gather(centerP/totalP);
end

plot(D, Pf);
[~, ii] = max(Pf);
fprintf('Optimal linear spacing: %.1f * MFD\n', spacings(ii));

return;
%% Test quadratic spacing
spacings = 2:0.1:6;
D = zeros(size(spacings));  Pf = D;

for ii = 1:numel(spacings)
    MFD = 5e-6;
    spacing = spacings(ii)*MFD;
    beams = arrayCPA(4);
    beams = beams.^2.*sign(beams)*0.5/max(beams(:))^2+beams;
    beams(:,1:2) = spacing*beams(:,1:2);
    x0 = 2*spacing*((size(beams,1)^0.5)/2+2); y0 = x0;

    [E0, x0, y0] = efieldGaussianBeam(x0, y0, beams, 'q', gaussianQ(lambda, MFD, z0), 'N', 2^9, 'plot', 0);
    totalP = sum((gradient(y0)*gradient(x0)) .* abs(E0).^2, 'all');

    [Ez, xz, yz] = efieldMeanKernel(x0, y0, z, E0, 'xz', xz, 'yz', yz, 'N', 2^8, 'valcheck', 0);

    centerR = 0.5/2;
    centerI = xz.^2 + yz.^2 < centerR^2;
    centerP = mean(abs(Ez(centerI)).^2, 'all')*pi*centerR^2;
    
    D(ii) = 2*max(beams(:,1));
    Pf(ii) = gather(centerP/totalP);
end

hold on; plot(D, Pf); hold off;
[~, ii] = max(Pf);
fprintf('Optimal quadratic spacing: (%.1f * MFD)^2\n', spacings(ii));

