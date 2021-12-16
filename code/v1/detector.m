function [F,T,P_time,P_freq] = detector(data,fs,Nfft,overlap,threshold,energy_perc)
%
% Author: Camille Pagniello (5/10/2019)
%
% INPUTS:
%   data - vector containing samples from .wav, int16
%   fs - sampling frequency of .wav
%   threshold - a scalar that you multiply by the mean noise level to get the threshold level
%   energy_perc - a scalar (0-100) that is used to set the threshold of the power filter

% OUTPUTS
%   P - all power spectra values
%   T - all times

%% Detection

%% Spectrogram
[~,F,T,P] = spectrogram(data,kaiser(Nfft,7.85),round(Nfft*overlap),Nfft,fs,'yaxis');
T = transpose(T);

%% Whitening
[Pxx,~] = pwelch(data,kaiser(Nfft,7.85),round(Nfft*overlap),Nfft,fs);
P = P./Pxx;

%% Power Threshold Filter
temp = P((((F < 400) | (F > 800))==0),:);
min_P = prctile(temp(:),energy_perc); % determine threshold for power filter
P(P < min_P) = 0; % apply dB filter

%% High Pass Filter
P(F < 50,:) = 0; % remove frequencies below 50 Hz

%% Band Pass Filter
P(F < 400,:) = 0;
P(F > 800,:) = 0;

%% Binarize Image, Fill Pieces and Remove Objects with Less than 25 Pixels
BW = imbinarize(10*log10(P));
BW = bwareaopen(BW,25);

P(BW == 0) = 0;

P_freq = P;

%% Time Series
P_time = transpose(sum(abs(P).^2));

% %% Noise Threshold
% avg_noise = mean(P(P >= 0)); % average noise level
% P(P < avg_noise*(threshold/100)) = 0;
% 
% %% Smooth
% P_time = smooth(P);

end