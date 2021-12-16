function [T,P] = camera_detector(data,fs,Nfft,overlap,freq,v,noise_threshold)
%
% Author: Camille Pagniello, cpagniel@ucsd.edu (last update: 1/29/2020)
%
% INPUTS:
%   data - vector containing samples from .wav, int16
%   fs - sampling frequency of .wav
%   v - power law exponents 
%   noise_threshold - a scalar (0-100) that is used to set the threshold of
%   the power filter to remove values below the noise level
%
% OUTPUTS
%   P - energy sum
%   T - all times

%% Spectrogram

[~,F,T,P] = spectrogram(data,kaiser(Nfft,7.85),round(Nfft*overlap),Nfft,fs,'yaxis');
T = transpose(T);

%% Frequency Standardization

P = P./(median(P,2));

%% Power Noise Threshold Filter

min_P = prctile(P(:),noise_threshold); % determine threshold for power filter
P(P < min_P) = 0; % apply dB filter

%% Limit Band to Selected Frequency Range

P(F < freq(1),:) = 0; % remove frequencies below 50 Hz
P(F > freq(2),:) = 0; % remove frequencies above 2000 Hz

%% Time Series

P = transpose(sum(P).^v); 

end