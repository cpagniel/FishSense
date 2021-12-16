%% Import Call

[data,fs] = audioread('test_07102018_230457_10s.wav','native');
overlap = 0.9; 
Nfft = 512;

%% Plot Original Call

figure()
% [~,F,T,P] = spectrogram(double(data(:,1))-mean(double(data(:,1))),kaiser(Nfft,7.85),round(Nfft*overlap),Nfft,fs,'yaxis');
imagesc(T,F,10*log10(P))
set(gca,'YDir','normal'); ylim([0 1500]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('Original Spectrogram');
colormap(jet); h = colorbar; caxis([0 15]); ylabel(h,{'Spectral Density','(dB re counts^2/Hz)'});

%% Whiten

% [Pxx,~] = pwelch(double(data(:,1))-mean(double(data(:,1))),kaiser(Nfft,7.85),round(Nfft*overlap),Nfft,fs);
% P = P./Pxx;

% P_Pxx(F < 500,:) = 0;
% P_Pxx(F > 1000,:) = 0;

figure()
imagesc(T,F,10*log10(P))
set(gca,'YDir','normal'); ylim([0 1500]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('Whitened Spectrogram');
colormap(jet); h = colorbar; caxis([0 15]); ylabel(h,{'Spectral Density','(dB re counts^2/Hz)'});

%% Threshold  

% threshold = prctile(P(:),97.5);
% P(P < threshold) = 0;

figure()
imagesc(T,F,10*log10(P))
set(gca,'YDir','normal'); ylim([0 1500]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('Whitened Spectrogram and 97.5% Threshold');
colormap(jet); h = colorbar; caxis([0 15]); ylabel(h,{'Spectral Density','(dB re counts^2/Hz)'});

%% Filtering

figure()
imagesc(T,F,10*log10(P))
set(gca,'YDir','normal'); ylim([0 1500]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('Whitened Spectrogram, 97.5% Threshold and Band Pass Filter 400-800 Hz');
colormap(jet); h = colorbar; caxis([0 15]); ylabel(h,{'Spectral Density','(dB re counts^2/Hz)'});

%% Power Spectrum

% P = sum(P,1);

figure()
plot(T,P);
% hold on
% plot(T,mean(P(P > 0))*0.70*ones(length(T),1))
