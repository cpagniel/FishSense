function TheFishSense(audiofile,out_dir,num_channels,sample_dur,call_duration,freq,v,noise_threshold,array)
%% The Fish Sense
% Generalized power-law detection algorithm for fish sounds based on Helble et al. (2012)
%
% INPUTS:
%   audiofile - string of complete path and name of .wav file (e.g., F:Detector\otherstuff\name.wav)
%   out_dir - string of output directory to save .txt file
%   channel - specify channel to use if multi-channel .wav file
%   sample_dur - duration of extracted sub sample from .wav file (seconds)
%   call_duration -
%   freq -
%   v - 
%   noise_threshold - a scalar (0-100) that is used to set the threshold of
%   the power filter to remove values below the noise level
% 
%
% Author: Camille Pagniello, cpagniel@ucsd.edu (last update: 1/29/2020)

%% Check Operating System

if ispc % PC
    slash='\';
elseif ismac % Mac
    slash='/';
else
    error('Invalid OS');
end

%% Initialize

disp('Initializing....');
tic; % start time for ETA

% pre-allocation
T = [];

%% Create Progress Bar

progress_bar = waitbar(0, 'Processing a new file...'); 

%% Analysis Parameters

fs = getfield(audioinfo(audiofile),'SampleRate'); % .wav sampling rate (Hz)
Nfft = 512; % FFT length
overlap = 0.9; % 90% overlap

%% Detection: Part 1

disp('Begin detection....');

start_sample = transpose(1:sample_dur*fs:getfield(audioinfo(audiofile),...
    'TotalSamples')); % start sample of each sub sample
end_sample = min(start_sample + sample_dur*fs - 1, getfield(audioinfo(audiofile),...
    'TotalSamples')); % end sample of each sub sample

for i = 1:(length(start_sample)-1)
    
    % Update Progress Bar
    elapsed_t = toc/60; % elapsed time (minutes)
    ETA = ((elapsed_t/(i))*(size(start_sample,1)-(i))); % ETA (minutes)
    waitbar(i/(size(start_sample,1)), progress_bar,... 
        sprintf('Current File Progress %d%%\n ETA: %d min',... 
        round(i/(size(start_sample,1))*100), round(ETA)));
    
    % Detector
    t = fish_detector(audiofile,start_sample(i),end_sample(i),num_channels,fs,Nfft,overlap,freq,v,noise_threshold,call_duration,array); 
    
    % Correct Times for Sub Sample Start
    for j = 1:length(t)
        t(j) = (start_sample(i)/fs)+t(j);
    end
    clear j
    
    % Concatenate
    T = vertcat(T,t);
    
    clear t
    clear data
end
clear i 

close(progress_bar);

%% Write to Table

var_names = {'StartSample','EndSample'};
out = table(T(:,1)*fs,T(:,2)*fs,'VariableNames',var_names);

[~,out_file,~]=fileparts(audiofile);
writetable(out,strcat(out_dir,slash,out_file,'_detect.csv'));

disp('Finish detection!');

toc

end


