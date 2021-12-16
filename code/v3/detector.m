function t = detector(app,audiofile,fs,start_sample,end_sample,hsens,tf)

%% Loop Through Channels
for i = 1:4 %app.NumberofChannels.Value
    
    %% Data
    [data,~] = audioread(audiofile,[start_sample end_sample],'native'); % load .wav file sub sample
    data = double(data(:,i));
    data = 10.^(-1*(hsens(i)+tf(i))/20)*data; % DAQ calibration
    data = data-mean(data); % remove mean
    
    % Spectrogram
    [~,F,T,P] = spectrogram(data,kaiser(app.FFTLength.Value,7.85),...
        round(app.FFTLength.Value*(app.PercentOverlap.Value/100)),app.FFTLength.Value,...
        fs,'yaxis');
    T = transpose(T);
    
    % Frequency and Time Standardization
    P = P./(median(P,2).*median(P,1));
    
    % Add Hydrophones Incoherently
    if i == 1
        P_total = zeros(size(P));
    end
    
    P_total = P_total + P;
    
end

P = P_total;

%% Power Noise Threshold Filter
min_P = prctile(P(:),app.NoiseThreshold.Value); % determine threshold for power filter
P(P < min_P) = 0; % apply dB filter

%% Limit Band to Selected Frequency Range
P(F < app.MinCallFrequency.Value,:) = 0; % remove frequencies below minimum frequency
P(F > app.MaxCallFrequency.Value,:) = 0; % remove frequencies above maximum frequency

%% Time Series
P = transpose(sum(P).^app.PowerLawExponent.Value);

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
    
    %% Discard Calls Less than Minimum Call Duration
    t(end_t-start_t < app.MinCallDuration.Value,:) = NaN;
    t(isnan(t(:,1)),:) = [];
    
     %% Discard Calls More than Maximum Call Duration
    t(end_t-start_t > app.MaxCallDuration.Value,:) = NaN;
    t(isnan(t(:,1)),:) = [];
    
else
    
    t = [];
    
end

end