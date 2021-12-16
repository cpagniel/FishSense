%% Plot Call

overlap = 0.9;
time = (0:length(data3)-1)./fs; % s
set(0,'DefaultFigureWindowStyle','normal')

figure('Position',[50 50 1850 1050])

subplot(441)
[~,F,T,P] = spectrogram(data1,kaiser(Nfft1,7.85),round(Nfft1*overlap),Nfft1,fs,'yaxis');
imagesc(T,F,10*log10(P))
set(gca,'YDir','normal'); ylim([yy(1) yy(2)]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title(['Spectrogram, Nfft = ' num2str(Nfft1)]);
colormap(jet); h = colorbar; caxis([-40 20]); ylabel(h,{'Spectrum level','(dB re counts^2/Hz)'});
subplot(442)
[~,F,T,P] = spectrogram(data1,kaiser(Nfft2,7.85),round(Nfft2*overlap),Nfft2,fs,'yaxis');
imagesc(T,F,10*log10(P))
set(gca,'YDir','normal'); ylim([yy(1) yy(2)]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title(['Spectrogram, Nfft = ' num2str(Nfft2)]);
colormap(jet); h = colorbar; caxis([-40 20]); ylabel(h,{'Spectrum level','(dB re counts^2/Hz)'});
subplot(443)
plot(time,filter_func(data1,fs,yy-1))
xlabel('Time (s)'); ylabel('Amplitude');
axis tight
subplot(444)
[Pxx,f] = pwelch(data1,kaiser(Nfft1,7.85),round(Nfft1*overlap),Nfft1,fs);
Pxx = 10*log10(Pxx);
plot(f,Pxx);
xlabel('Frequency (Hz)'); ylabel({'Power spectral density';'(dB re counts^2/Hz)'}); title(['Spectrum, Nfft = ' num2str(Nfft1)]);
xlim([yy(1) yy(2)]);

subplot(445)
[~,F,T,P] = spectrogram(data2,kaiser(Nfft1,7.85),round(Nfft1*overlap),Nfft1,fs,'yaxis');
imagesc(T,F,10*log10(P))
set(gca,'YDir','normal'); ylim([yy(1) yy(2)]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title(['Spectrogram, Nfft = ' num2str(Nfft1)]);
colormap(jet); h = colorbar; caxis([-40 20]); ylabel(h,{'Spectrum level','(dB re counts^2/Hz)'});
subplot(446)
[~,F,T,P] = spectrogram(data2,kaiser(Nfft2,7.85),round(Nfft2*overlap),Nfft2,fs,'yaxis');
imagesc(T,F,10*log10(P))
set(gca,'YDir','normal'); ylim([yy(1) yy(2)]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title(['Spectrogram, Nfft = ' num2str(Nfft2)]);
colormap(jet); h = colorbar; caxis([-40 20]); ylabel(h,{'Spectrum level','(dB re counts^2/Hz)'});
subplot(447)
plot(time,filter_func(data2,fs,yy-1))
xlabel('Time (s)'); ylabel('Amplitude');
axis tight
subplot(448)
[Pxx,f] = pwelch(data2,kaiser(Nfft1,7.85),round(Nfft1*overlap),Nfft1,fs);
Pxx = 10*log10(Pxx);
plot(f,Pxx);
xlabel('Frequency (Hz)'); ylabel({'Power spectral density';'(dB re counts^2/Hz)'}); title(['Spectrum, Nfft = ' num2str(Nfft1)]);
xlim([yy(1) yy(2)]);

subplot(449)
[~,F,T,P] = spectrogram(data3,kaiser(Nfft1,7.85),round(Nfft1*overlap),Nfft1,fs,'yaxis');
imagesc(T,F,10*log10(P))
set(gca,'YDir','normal'); ylim([yy(1) yy(2)]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title(['Spectrogram, Nfft = ' num2str(Nfft1)]);
colormap(jet); h = colorbar; caxis([-40 20]); ylabel(h,{'Spectrum level','(dB re counts^2/Hz)'});
subplot(4,4,10)
[~,F,T,P] = spectrogram(data3,kaiser(Nfft2,7.85),round(Nfft2*overlap),Nfft2,fs,'yaxis');
imagesc(T,F,10*log10(P))
set(gca,'YDir','normal'); ylim([yy(1) yy(2)]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title(['Spectrogram, Nfft = ' num2str(Nfft2)]);
colormap(jet); h = colorbar; caxis([-40 20]); ylabel(h,{'Spectrum level','(dB re counts^2/Hz)'});
subplot(4,4,11)
plot(time,filter_func(data3,fs,yy-1))
xlabel('Time (s)'); ylabel('Amplitude');
axis tight
subplot(4,4,12)
[Pxx,f] = pwelch(data3,kaiser(Nfft1,7.85),round(Nfft1*overlap),Nfft1,fs);
Pxx = 10*log10(Pxx);
plot(f,Pxx);
xlabel('Frequency (Hz)'); ylabel({'Power spectral density';'(dB re counts^2/Hz)'}); title(['Spectrum, Nfft = ' num2str(Nfft1)]);
xlim([yy(1) yy(2)]);

subplot(4,4,13)
[~,F,T,P] = spectrogram(data4,kaiser(Nfft1,7.85),round(Nfft1*overlap),Nfft1,fs,'yaxis');
imagesc(T,F,10*log10(P))
set(gca,'YDir','normal'); ylim([yy(1) yy(2)]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title(['Spectrogram, Nfft = ' num2str(Nfft1)]);
colormap(jet); h = colorbar; caxis([-40 20]); ylabel(h,{'Spectrum level','(dB re counts^2/Hz)'});
subplot(4,4,14)
[~,F,T,P] = spectrogram(data4,kaiser(Nfft2,7.85),round(Nfft2*overlap),Nfft2,fs,'yaxis');
imagesc(T,F,10*log10(P))
set(gca,'YDir','normal'); ylim([yy(1) yy(2)]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title(['Spectrogram, Nfft = ' num2str(Nfft2)]);
colormap(jet); h = colorbar; caxis([-40 20]); ylabel(h,{'Spectrum level','(dB re counts^2/Hz)'});
subplot(4,4,15)
plot(time,filter_func(data4,fs,yy-1))
xlabel('Time (s)'); ylabel('Amplitude');
axis tight
subplot(4,4,16)
[Pxx,f] = pwelch(data4,kaiser(Nfft1,7.85),round(Nfft1*overlap),Nfft1,fs);
Pxx = 10*log10(Pxx);
plot(f,Pxx);
xlabel('Frequency (Hz)'); ylabel({'Power spectral density';'(dB re counts^2/Hz)'}); title(['Spectrum, Nfft = ' num2str(Nfft1)]);
xlim([yy(1) yy(2)]);

cd(Model.resultfn)
mkdir([labels{jj} '_' datestr(cname,'yymmddHHMMSSFFF') '.wav'])
cd([labels{jj} '_' datestr(cname,'yymmddHHMMSSFFF') '.wav'])

saveas(gcf,'call.fig')
saveas(gcf,'call.png')

disp('data plotted')