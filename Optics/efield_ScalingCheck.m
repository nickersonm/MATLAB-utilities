lambda = 1.55e-6; k = 2*pi/lambda;
MFD = 5e-6;
sz = MFD*4;

%%
powG = [];
Ns = 4:12;
for ii=Ns
    [E, x, y] = efieldGaussianBeam(sz, sz, [0 0], 'q', gaussianQ(lambda, MFD, 0), 'N', 2^ii, 'plot', 0);
%     gridscale = sqrt(gradient(y) .* gradient(x));
    powG(end+1) = gather(mean(abs(E).^2,'all')*sz^2) *pi/(MFD^2);
%     powG(end+1) = gather(sum(abs(gridscale.*E).^2,'all')) /(MFD^2*pi/2);
end

figure(1);
plot(Ns,powG);


%%
powMK = [];
Ns = 4:12;
E0 = efieldGaussianBeam(sz, sz, [0 0], 'q', gaussianQ(lambda, MFD, 0), 'N', 2^9, 'plot', 0);
P0 = gather(mean(abs(E0).^2,'all')*sz^2);
for ii=Ns
    E = efieldMeanKernel(sz, sz, 1e-5, E0, 'plot', 0, 'N', 2^ii);
    powMK(end+1) = gather(mean(abs(E).^2,'all')*sz^2)/P0;
end

figure(1);
hold on; plot(Ns, powMK); hold off;


%%
powMK = [];
scales = 0:0.5:4; scales = 2*sqrt(scales)+1;
E0 = efieldGaussianBeam(sz, sz, [0 0], 'q', gaussianQ(lambda, MFD, 0), 'N', 2^9, 'plot', 0);
P0 = gather(mean(abs(E0).^2,'all')*sz^2);
for sc=scales
    E = efieldMeanKernel(sz, sz, sz, E0, 'plot', 0, 'N', 2^8, 'xz', sz*sc, 'yz', sz*sc);
    powMK(end+1) = gather(mean(abs(E).^2,'all')*(sz*sc)^2)/P0;
end

figure(2);
plot(scales, powMK);

%%
E0 = efieldGaussianBeam(sz, sz, [0 0], 'q', gaussianQ(lambda, MFD, 0), 'N', 2^9);
efieldMeanKernel(sz, sz, 100*MFD, E0, 'plot', 3, 'N', 2^8, 'xz', 200*MFD, 'yz', 200*MFD, 'valcheck', 1);
