% Filter testing
%   Michael Nickerson 2025-03-21


%% Defaults
normal = true;
dt = 1e-10; 
analogFreq = 1/dt/5;
analogOrder = 3;
Nfft = 16;
F = 1/dt;


%% Set up time series
t = dt*2^16;
t = linspace(0, t, 2^16);
x0 = randn(size(t));


%% Generate and apply filters
% Generate analog filter function
% Assemble filter: low-pass butterworth
afilt = designfilt("lowpassiir", "DesignMethod", "butter", ...
                   "FilterOrder", analogOrder, ...
                   "HalfPowerFrequency", 2*analogFreq/F);

[b, a] = butter(analogOrder, 2*pi*analogFreq/F, 's');
[ba, aa] = bilinear(b, a, 1);

x1 = filtfilt(afilt, x0);
x2 = filter(ba, aa, x0);


%% FT
fft0 = fftshift(abs(fft(x0, Nfft*numel(x0))).^2); fft0 = fft0(end/2:end);
fft1 = fftshift(abs(fft(x1, Nfft*numel(x1))).^2); fft1 = fft1(end/2:end);
fft2 = fftshift(abs(fft(x2, Nfft*numel(x2))).^2); fft2 = fft2(end/2:end);


%% Plot
f = linspace(0, 0.5, numel(fft0))/dt;

figure(1);
plot(f/1e9, fft0, f/1e9, fft1, f/1e9, fft2);
xlim([min(f) 2*analogFreq]/1e9); grid on;
legend(["Original", "`designfilt`", "`butter`"], "Location", "ne", "FontSize", 14);
xlabel("Frequency [GHz]", "FontWeight","bold");
ylabel("Amplitude [arb]", "FontWeight","bold")



%% Arbitrary frequency-domain shape
f3 = (square(t*9*pi/max(t))+1)/2;
x3 = ifft(fft(x0) .* f3);

fft3 = fftshift(abs(fft(x3, Nfft*numel(x3))).^2); fft3 = fft3(end/2:end);

% Plot
figure(2);
plot(f/1e9, fft0, f/1e9, fft3, linspace(0, 0.5, numel(x0)/2)/dt/1e9, f3(1:end/2)*max(fft0)/max(f3));
xlim([min(f) 2*analogFreq]/1e9); grid on;
legend(["Original", "`ifft`", "Expected Response"], "Location", "ne", "FontSize", 14);
xlabel("Frequency [GHz]", "FontWeight","bold");
ylabel("Amplitude [arb]", "FontWeight","bold")
