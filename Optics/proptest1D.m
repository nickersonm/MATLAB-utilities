%% proptest1D.m
%   Compare H-F and Fraunhofer 1D propagation

clear; %close("all");


%% Assemble a near-field pattern
lambda = 1.55e-6;
N = 2^12;

nTx     = 128;
dx      = 1.5e-6;
x = (1:nTx) * dx; x = x(:) - mean(x);

steer = 10*pi/180;
th = steer + [-1 1] * 8*lambda/nTx/dx; th = linspace(min(th), max(th), N);
ph = 2*pi*x/lambda * sin(steer);
E0 = ones(size(x)) .* exp(1i*ph) / numel(x).^2; % Uniform emission


%% Propagate to the far field
[EzHF, thHF] = simpleHuygensFresnel1D(x, E0, "lambda", lambda, "th", th);
[EzFF, thFF] = simpleFraunhofer1D(x, E0, "lambda", lambda, "th", th);


%% Reverse propagation
[E0HF, xHF] = simpleHuygensFresnel1D(thHF, EzHF, "lambda", lambda, "th", x, "reverse");
[E0FF, xFF] = simpleFraunhofer1D(thHF, EzHF, "lambda", lambda, "th", x, "reverse");


%% Plot
figureSize(1, 1400, 600); clf;
mgn = [0.12, 0.05];
h = subplot_tight(1,2,1, mgn);
set(gca, "Position", gca().Position - [0.01 0 0 0]);
plot(thFF*180/pi, 10*log10(abs(EzFF).^2), "LineWidth", 2); hold on;
plot(thHF*180/pi, 10*log10(abs(EzHF).^2), "--", "LineWidth", 2);
ll = ["Fraunhofer FFT", "Huygens-Fresnel"];
hold off; grid on; axis tight;
xlabel("Angle [°]");
ylabel("Intensity [dB]");
legend(ll, "Location", "nw");
title(h, "Forward", "FontSize", 14);

h = subplot_tight(1,2,2, mgn);
set(gca, "Position", gca().Position + [0.01 0 0 0]);
plot(xFF, 10*log10(abs(E0FF).^2), "LineWidth", 2); hold on;
plot(xHF, 10*log10(abs(E0HF).^2), "--", "LineWidth", 2);
plot(x, ones(size(E0)) * mean(10*log10(abs(E0HF).^2)), "LineWidth", 2);
ll = ["Fraunhofer FFT", "Huygens-Fresnel", "Original (not scaled)"];
hold off; grid on; axis tight;
xlabel("Position [m]");
ylabel("Intensity [dB]");
legend(ll, "Location", "s");
title(h, "Backward", "FontSize", 14);

% figure(2); clf;
% plot(xFF, angle(E0FF), "LineWidth", 2); hold on;
% plot(xHF, angle(E0HF), "LineWidth", 2); hold off;
