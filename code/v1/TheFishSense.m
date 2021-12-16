function TheFishSense(audiofile, out_dir, channel, sub_sample_dur, threshold, energy_perc)
%% The Fish Sense
% Detector for loud sounds. Optimized for UF1 kelp forest calls.
%
% INPUTS:
%   audiofile - string of complete path and name of .wav file (e.g., F:Detector\otherstuff\name.wav)
%   out_dir - string of output directory to save .txt file
%   channel - specify channel to use if multi-channel .wav file
%   sub_sample_dur - duration of extracted sub sample from .wav file (seconds)
%   threshold - a scalar that you multiply by the mean noise level to get the threshold level
%   energy_perc - a scalar (0-100) that is used to set the threshold of the power filter
%
% Author: Camille Pagniello (5/10/2019)

%% Operating System

if ispc
    slash='\';
elseif ismac
    slash='/';
else
    error('Invalid OS');
end

%% Initialize

disp('Initialize....');
tic; % start time for ETA

% Pre-Allocation
P_t = [];
P_f = [];
T = [];

% Create Progress Bar
progress_bar = waitbar(0, 'Processing a new file...'); 

% Analysis Parameters
fs = getfield(audioinfo(audiofile),'SampleRate'); % .wav sampling rate (Hz)
Nfft = 512; % FFT length
overlap = 0.9; % 90% overlap

%% Detection

disp('Begin detection....');

%% Run Detector on Entire .wav File

start_sample = transpose(1:sub_sample_dur*fs:getfield(audioinfo(audiofile),...
    'TotalSamples')); % start sample of each sub sample
end_sample = min(start_sample + sub_sample_dur*fs - 1, getfield(audioinfo(audiofile),...
    'TotalSamples')); % end sample of each sub sample

for i = 1:(length(start_sample)-1)
    [data,~] = audioread(audiofile, [start_sample(i) end_sample(i)],'native'); % load .wav file sub sample
    data = double(data(:,channel));
    data = 10.^(78.2/20)*data; % DAQ calibration
    data = data-mean(data); % remove mean
    
    % Update Progress Bar
    elapsed_t = toc/60; % elapsed time (minutes)
    ETA = ((elapsed_t/(i))*(size(start_sample,1)-(i))); % ETA (minutes)
    waitbar(i/(size(start_sample,1)), progress_bar,... 
        sprintf('Current File Progress %d%%\n ETA: %d min',... 
        round(i/(size(start_sample,1))*100), round(ETA)));
    
    % Detection
    [f,t,p_t,p_f] = detector(data, fs, Nfft, overlap, threshold, energy_perc); 
    
    % Correct Times for Sub Sample Start
    for j = 1:length(t)
        t(j) = (start_sample(i)/fs)+t(j);
    end
    
    % Concatenate
    P_t = vertcat(P_t,p_t);
    P_f = horzcat(P_f,p_f);
    T = vertcat(T,t);
    
end

F = f;

clear i j p_t p_f f t

close(progress_bar);

%% Extract Calls

% Noise Threshold
avg_noise = mean(P_t(P_t >= 0)); % average noise level
P_t(P_t < avg_noise*(threshold/100)) = 0;

% Smooth
P_t = smooth(P_t);

% Start and End Times
start_t = T(diff(P_t ~= 0)== 1);
end_t = T(find(diff(P_t ~= 0)== -1)+1);
if length(start_t) > length(end_t)
    start_t(end) = [];
elseif length(end_t) > length(start_t)
    end_t(1) = [];
end
t = horzcat(start_t,end_t);

% Discard calls less than 0.25 seconds and more than 0.70 seconds
dur = end_t-start_t;
t(dur < 0.25,:) = NaN; t(dur > 0.70,:) = NaN;
t(isnan(t(:,1)),:) = [];

%% Write to Table

var_names = {'StartSample','EndSample'};
out = table(t(:,1)*fs,t(:,2)*fs,'VariableNames',var_names);

[~,out_file,~]=fileparts(audiofile);
writetable(out,strcat(out_dir,slash,out_file,'_detect.csv'));

disp('Finish detection!');

toc

end


