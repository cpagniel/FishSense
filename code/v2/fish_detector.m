function t = fish_detector(audiofile,start_sample,end_sample,num_channels,fs,Nfft,overlap,freq,v,noise_threshold,call_duration,array)
%
% Author: Camille Pagniello, cpagniel@ucsd.edu (last update: 1/29/2020)
%
% INPUTS:
%   data - vector containing samples from .wav, int16
%   fs - sampling frequency of .wav
%   v - power law exponents
%   noise_threshold - a scalar (0-100) that is used to set the threshold of
%   the power filter to remove values below the noise level
%   call_duration - minimum call duration (seconds)
%
% OUTPUTS
%   P - energy sum
%   T - all times

%% Loop Through Channels

for i = 1:num_channels
    
    %% Data
    
    [data,~] = audioread(audiofile, [start_sample end_sample],'native'); % load .wav file sub sample
    data = double(data(:,i));
    data = 10.^(-1*(array.hsens(i)+array.tf)/20)*data; % DAQ calibration
    data = data-mean(data); % remove mean

    %% Spectrogram

    [~,F,T,P] = spectrogram(data,kaiser(Nfft,7.85),round(Nfft*overlap),Nfft,fs,'yaxis');
    T = transpose(T);

    %% Frequency and Time Standardization

    P = P./(median(P,2).*median(P,1));
    
    %% Add hydrophones (incoherent noise)
    
    if i == 1
        P_total = zeros(size(P));
    end
    
    P_total = P_total + P;
    
end

P = P_total;

%% Power Noise Threshold Filter

min_P = prctile(P(:),noise_threshold); % determine threshold for power filter
P(P < min_P) = 0; % apply dB filter

%% Limit Band to Selected Frequency Range

P(F < freq(1),:) = 0; % remove frequencies below 50 Hz
P(F > freq(2),:) = 0; % remove frequencies above 1000 Hz

%% Time Series

P = transpose(sum(P).^v);

%% Detection Threshold

detection_threshold = mean(P); % mean of incoherent sum over frequency of spectrogram

%% Start and End Indices

start_ind = find(diff(P ~= 0)== 1);
end_ind = find(diff(P ~= 0)== -1);

if ~isempty(start_ind) && ~isempty(end_ind)
    if end_ind(1) < start_ind(1)
        start_ind = [1; start_ind];
    end

    if end_ind(end) < start_ind(end)
        end_ind = [end_ind; length(P)];
    end

    if length(start_ind) > length(end_ind)
        start_ind(end) = [];
    elseif length(end_ind) > length(start_ind)
        end_ind(1) = [];
    end

%% Signal Present or Absent

keep = [];
for i = 1:length(start_ind)
    d = sum(P(start_ind(i):end_ind(i)) > detection_threshold);
    if d > 0
        keep = [keep; i];
    end
end

%% Start and End Times

start_t = T(start_ind(keep));
end_t = T(end_ind(keep));

t = horzcat(start_t,end_t);

%% Discard calls less than call_duration

t(end_t-start_t < call_duration,:) = NaN;
t(isnan(t(:,1)),:) = [];

else
    
    t = [];

end

end