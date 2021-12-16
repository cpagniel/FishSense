%% Plotting

%% Spectrogram

figure()
imagesc(T,F,10*log10(P))
set(gca,'YDir','normal'); xlim([0 max(T)]); ylim([0 1500]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); 
colormap(jet); h = colorbar; ylabel(h,{'Spectral Density (dB re \muPa^2/Hz)'});
caxis([-10 10]);

%% Whitening

imagesc(T,F,10*log10(P))
set(gca,'YDir','normal'); xlim([0 max(T)]); ylim([0 1500]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('Whitened Spectrogram');
colormap(jet); h = colorbar; caxis([0 15]); ylabel(h,{'Spectral Density (dB re counts^2/Hz)'});

%% High Pass Filter and dB or Power Threshold Filter

imagesc(T,F,10*log10(P))
set(gca,'YDir','normal'); xlim([0 max(T)]); ylim([0 1500]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('Whitened Spectrogram and 97.5% Threshold');
colormap(jet); h = colorbar; caxis([0 15]); ylabel(h,{'Spectral Density(dB re counts^2/Hz)'});

%% Band Pass Filter and Binarize Image and Remove Objects with Less than 15 Pixels

imagesc(T,F,10*log10(P))
set(gca,'YDir','normal'); xlim([0 max(T)]); ylim([0 1500]);
xlabel('Time (s)'); ylabel('Frequency (Hz)');  title('Whitened Spectrogram, 97.5% Threshold, Band Pass Filter 400-800 Hz and Objects Less Than 50 Pixels Removed');
colormap(jet); h = colorbar; caxis([0 15]); ylabel(h,{'Spectral Density (dB re counts^2/Hz)'});

%% Power Spectra

figure()
plot(T,P_time)
xlabel('Time (s)'); ylabel('sum(|P|^2) (dB)'); %title('Power Spectra');

%% Noise Threshold and Smooth

plot(T,P)
xlabel('Time (s)'); ylabel('Power (dB)'); title('Power Spectra with Noise Threshold and Smoothed');
