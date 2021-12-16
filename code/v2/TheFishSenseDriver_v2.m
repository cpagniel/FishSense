%% The Fish Driver
% The Fish Driver runs: 
%   1 - TheFishSense: generalized power-law detection algorithm for fish
%   sounds on all channels
%   2 - TheCameraSense: generalized power-law detection algorithm for
%   camera shutter on all channels
%   3 - TheCameraCheck: checks to see if fish detections occur within 10
%   seconds of camera shutter
%   4 - TheFishFinder: localization of fish sound
%   5 -  

% User must define detector parameters before running and specify location
% of file containing array geometry, hydrophone sensitivity and camera
% start times.

% Author: Camille Pagniello, cpagniel@ucsd.edu (last update: 2/1/2020)

%% Define Detector Parameters

array.file = 'F:\C_files\C_Winter_2020\july2018_hydrophone_locations.txt'; % file containing array geometry and hydrophone sensitivity

% For Fish Calls
fish.sample_dur = 10; % duration of extracted sub sample from .wav file (seconds)
fish.noise_threshold = 95; % a scalar (0-100) that is used to set the threshold of the power filter to remove values below the noise level
fish.num_channels = 4; % number of channels to process
fish.freq = [50 1000]; % frequency (Hz)
fish.v = 2.5; % power law exponent, v = 2.5 is the optimal value in many cases from Nuttal's results
fish.call_duration = 0.1; % minimum duration of calls (seconds)

% For Camera Shutter
camera.sample_dur = 10; % duration of extracted sub sample from .wav file (seconds)
camera.noise_threshold = 99; % a scalar (0-100) that is used to set the threshold of the power filter to remove values below the noise level
camera.num_channels = 4; % number of channels to process
camera.freq = [50 2000]; % frequency (Hz)
camera.v = 2.5; % power law exponent, v = 2.5 is the optimal value in many cases from Nuttal's results
camera.call_duration = [0.24 0.27]; % [min max] minimum and maximum duration of shutter for channels 1, 2 and 4 (seconds)
camera.cam2ch = 3; % camera shutter is different on camera 2

%% Check Operating System

if ispc % PC
    slash='\';
elseif ismac % Mac
    slash='/';
else
    error('Invalid OS');
end

%% Array Geometry

array.M = 4; % number of hydrophones

array.xr = zeros(array.M,1); array.yr = zeros(array.M,1); array.zr = zeros(array.M,1); % hydrophone positions in x, y, and z
array_raw = load(array.file);
array.GPS = array_raw(:,1:3);

for m = 1:array.M
    [array.yr(m),array.xr(m)] = latlon2xy(array.GPS(m,1),array.GPS(m,2),...
        array.GPS(1,1),array.GPS(1,2));
    array.zr(m) = array.GPS(m,3);
end
array.xr = -1*array.xr;
array.zr = distdim(array.zr,'ft','m');
clear m

array.hsens = array_raw(:,4); % hydrophone sensitivity for channels 1 to 4 (from calibration sheet)
array.tf = 20*log10(2^16/(1.5--1.5)); % transfer function for ST4300 DAQ

for i = 2:array.M
    array.dist(i-1) = sqrt((array.xr(1)-array.xr(i))^2 + (array.yr(1)-array.yr(i))^2 + (array.zr(1)-array.zr(i))^2);
    array.tmax(i-1) = array.dist(i-1)/1500;
end

%% Camera Start Times

camera.time(1) = datenum(array_raw(1,7),array_raw(1,5),array_raw(1,6),array_raw(1,8),array_raw(1,9),array_raw(1,10)); 
camera.time(2) = datenum(array_raw(2,7),array_raw(2,5),array_raw(2,6),array_raw(2,8),array_raw(2,9),array_raw(2,10)); 
camera.time(3)  = datenum(array_raw(3,7),array_raw(3,5),array_raw(3,6),array_raw(3,8),array_raw(3,9),array_raw(3,10)); 
camera.time(4) = datenum(array_raw(4,7),array_raw(4,5),array_raw(4,6),array_raw(4,8),array_raw(4,9),array_raw(4,10)); 

clear array_raw

%% Get .wav Input File Directories

in_dir = uigetdir(pwd,'Select directory of .wav files');
files = dir(strcat(in_dir,slash,'*.wav')); % structure that contains information about the files to be processed

% Check for Additional Directories
choice = 'YES';
while strcmp(choice,'YES')
    choice = questdlg('Do you want to add additional directories to process?','Select Another Directory','YES','NO','YES');
    if strcmp(choice,'YES')
        in_dir = uigetdir(pwd,'Select additional directory of .wav files');
        files = [files; dir(strcat(in_dir,slash,'*.wav'))];
    end
end

%% Define Output File Directories

fish.out_dir = uigetdir(pwd,'Select directory to save output .csv files for initial fish call detections');
camera.out_dir = uigetdir(pwd,'Select directory to save output .csv files for camera shutter detections');
fish.out_dir2 = uigetdir(pwd,'Select directory to save output .csv files for fish calls detected on all four channels that are near a camera shutter');

%% TheFishSense

% Create Progress Bar for Fish Calls Initial Detection
dur = zeros(size(files,1),1); % pre-allocating
for i = 1:length(files)
    dur(i) = getfield(audioinfo(fullfile(files(i).folder,files(i).name)),'Duration'); % file duration (seconds)
end
disp(['Total Duration of Data = ',num2str(round(sum(dur)/(60*60*24))),' days']);

progress_bar_fish = waitbar(0, 'Processing data set...','Position',[500 500 275 60]);
proc_start_t = tic;

% Run Detector for Fish Calls
for i = 1:length(files)
    if datenum(files(i).name(11:end-8),'yymmddHHMMSS') > min(camera.time)
        if str2double(files(i).name(17:18)) >= 4 && str2double(files(i).name(17:18)) <= 22
            TheFishSense(fullfile(files(i).folder,files(i).name),fish.out_dir,...
                fish.num_channels,fish.sample_dur,fish.call_duration,fish.freq,...
                fish.v,fish.noise_threshold,array);
        end
    end
    
    % Update Progress Bar
    elapsed_t = toc(proc_start_t)/(60); % elapsed time (minutes)
    ETA = ((elapsed_t/sum(dur(1:i)))*(sum(dur)-sum(dur(1:i))));
    waitbar((sum(dur(1:i))/sum(dur)),...
        progress_bar_fish,sprintf('Overall Progress: %d%%\n ETA: %d min',...
        round((sum(dur(1:i))/sum(dur))*100),round(ETA)));
end

close(progress_bar_fish);
disp('Done Detection!');

%% TheCameraSense

% Create Progress Bar for Camera Shutter
dur = zeros(size(files,1),1); % pre-allocating
for i = 1:length(files)
    dur(i) = getfield(audioinfo(fullfile(files(i).folder,files(i).name)),'Duration'); % file duration (seconds)
end
disp(['Total Duration of Data = ',num2str(round(sum(dur)/(60*60*24))),' days']);

progress_bar_camera = waitbar(0, 'Finding camera shutters...','Position',[500 500 275 60]);
proc_start_t = tic;

% Run Detector for Camera Shutter
for i = 1:length(files)
    for channel = 1:camera.num_channels %getfield(audioinfo(fullfile(files(i).folder,files(i).name)),'NumChannels')
        if datenum(files(i).name(11:end-8),'yymmddHHMMSS') > camera.time(channel)
            if str2double(files(i).name(17:18)) >= 4 && str2double(files(i).name(17:18)) <= 22
                TheCameraSense(fullfile(files(i).folder,files(i).name),camera.out_dir,channel,...
                    camera.sample_dur,camera.call_duration,camera.freq,...
                    camera.v,camera.noise_threshold,camera.cam2ch,array);
            end
        end
    end
    
    % Update Progress Bar
    elapsed_t = toc(proc_start_t)/(60); % elapsed time (minutes)
    ETA = ((elapsed_t/sum(dur(1:i)))*(sum(dur)-sum(dur(1:i))));
    waitbar((sum(dur(1:i))/sum(dur)),...
        progress_bar_camera,sprintf('Overall Progress: %d%%\n ETA: %d min',...
        round((sum(dur(1:i))/sum(dur))*100),round(ETA)));
end

close(progress_bar_camera);
disp('Done Detection Camera Shutter!');

%% TheCameraCheck

% Create Progress Bar for CameraCheck
dur = zeros(size(files,1),1); % pre-allocating
for i = 1:length(files)
    dur(i) = getfield(audioinfo(fullfile(files(i).folder,files(i).name)),'Duration'); % file duration (seconds)
end
disp(['Total Duration of Data = ',num2str(round(sum(dur)/(60*60*24))),' days']);

progress_bar_camera = waitbar(0, 'Checking if detections are on all channels...','Position',[500 500 275 60]);
proc_start_t = tic;

% Check If Detections Are Within +/- 5 seconds of a Camera Shutter
for i = 1:length(files)
    if datenum(files(i).name(11:end-8),'yymmddHHMMSS') > min(camera.time)
            if str2double(files(i).name(17:18)) >= 4 && str2double(files(i).name(17:18)) <= 22
                TheCameraCheck(fullfile(files(i).folder,files(i).name),fullfile(camera.out_dir,files(i).name),fullfile(fish.out_dir,files(i).name),fish.out_dir2,array,camera);
            end
    end
    
        % Update Progress Bar
    elapsed_t = toc(proc_start_t)/(60); % elapsed time (minutes)
    ETA = ((elapsed_t/sum(dur(1:i)))*(sum(dur)-sum(dur(1:i))));
    waitbar((sum(dur(1:i))/sum(dur)),...
        progress_bar_camera,sprintf('Overall Progress: %d%%\n ETA: %d min',...
        round((sum(dur(1:i))/sum(dur))*100),round(ETA)));
end

close(progress_bar_camera);

%%

clear
