% audio signal processing
% question.2
clc;
[audio, Fs] = audioread('bird_sound.mp3');

% plot the signal in time domain
plot([1:length(audio)]./Fs,audio);
xlim([0 length(audio)/Fs]);
title('Song Signal','interpreter','latex');
xlabel('Time(s)','interpreter','latex');
ylabel('Amplitude','interpreter','latex');
grid on;

% power spectrum of the signal
figure
X = audio;
L = length(X);
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
plot(f,P1);
xlim([0 max(f)])
title('Fast Fourier Transform Of The AM Signal','interpreter','latex');
xlabel('Frequency(Hz)','interpreter','latex');
ylabel('Amplitude','interpreter','latex');
grid on;

% if we look at the power spectrum, we relize we have the most power in the
% frequency range of [2k, 7k] Hz

% calculate the power in this range using bandpower 
SelectedRangePower = bandpower(audio,Fs,[2000 7000]);

% calculate the total power of the signal
TotalPower = bandpower(audio,Fs,[0 max(f)]);

% now if calculate ratio of a to b, it will be a big percentage showing
% that we have most of the signal`s energy in the selected frequency range
Percentage = (SelectedRangePower/TotalPower)*100;


% remove AM modulation and export the bird sound
[yupper,ylower] = envelope(audio);
% power spectrum of the madulated signal
figure
X = ylower+mean(ylower);
L = length(X);
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
plot(f,P1);

title('Fast Fourier Transform Of The Message Signal','interpreter','latex');
xlabel('Frequency(Hz)','interpreter','latex');
ylabel('Amplitude','interpreter','latex');
grid on;