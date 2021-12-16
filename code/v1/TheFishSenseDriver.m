%% The Fish Sense Driver
% This driver can process multi-channel .wav files from multiple folders, 
% across multiple locations. It calls TheFishSense function that is tuned
% to detect UF1 type calls from the kelp forest. Currently hard set to only
% process the first channel.

% A .csv file is saved in a specified location for each .wav file with
% detections.
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

%% File Directory

f_dir = uigetdir(pwd,'Select directory of .wav files');
files = dir(strcat(f_dir,slash,'*.wav')); % structure that contains information about the files to be processed

% Additional Directories
choice='YES';
while strcmp(choice,'YES')
    choice = questdlg('Do you want to add additional directories to process?','Select Another Directory','YES','NO','YES');
    if strcmp(choice,'YES')
        f_dir = uigetdir(pwd,'Select additional directory of .wav files');
        files = [files; dir(strcat(f_dir,slash,'*.wav'))];
    end
end

%% Output Directory

out_dir = uigetdir(pwd,'Select directory to save output .txt files');

%% Progress Bar

dur = zeros(size(files,1),1); % pre-allocating
for i = 1:length(files)
    dur(i) = getfield(audioinfo(fullfile(files(i).folder,files(i).name)),'Duration'); % file duration (seconds)
end
disp(['Total Duration = ',num2str(round(sum(dur)/(60*60*24))),' days']);

%% Run Detector

% Progress Bar
progress_bar_main = waitbar(0, 'Processing Data Set...','Position',[500 500 275 60]);
proc_start_t = tic;

% Detector Parameters
sub_sample_dur = 5; % seconds
threshold = 100; % above average noise
energy_perc = 90; % percent
channel = 1;

for i = 1:length(files)
    % for channel = 1:getfield(audioinfo(fullfile(files(i).folder,files(i).name)),'NumChannels')
        TheFishSense(fullfile(files(i).folder,files(i).name),out_dir,channel,sub_sample_dur,threshold,energy_perc);
    % end
    
    % Progress Bar
    elapsed_t = toc(proc_start_t)/(60); % elapsed time (minutes)
    ETA = ((elapsed_t/sum(dur(1:i)))*(sum(dur)-sum(dur(1:i))));
    waitbar((sum(dur(1:i))/sum(dur)), progress_bar, ...
        sprintf('Overall Progress: %d%%\n ETA: %d min',...
        round((sum(dur(1:i))/sum(dur))*100),round(ETA)));
end

close(progressBar_main);
disp('Done!');

clear
