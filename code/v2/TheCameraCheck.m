function TheCameraCheck(audiofile,camerafile,detectfile,out_dir,array,camera)
%% The Camera Check
% 
%
% INPUTS:

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

%% Load .csv Camera File

if datenum(camerafile(61:end-8),'yymmddHHMMSS') > min(camera.time(1))
    if exist([camerafile(1:end-4) '_ch1_detect.csv'],'file')
        cam_1 = table2array(readtable([camerafile(1:end-4) '_ch1_detect.csv']),'Delimiter',',');
    end
end
if datenum(camerafile(61:end-8),'yymmddHHMMSS') > min(camera.time(2))
    if exist([camerafile(1:end-4) '_ch2_detect.csv'],'file')
        cam_2 = table2array(readtable([camerafile(1:end-4) '_ch2_detect.csv']),'Delimiter',',');
    end
end
if datenum(camerafile(61:end-8),'yymmddHHMMSS') > min(camera.time(3))
    if exist([camerafile(1:end-4) '_ch3_detect.csv'],'file')
        cam_3 = table2array(readtable([camerafile(1:end-4) '_ch3_detect.csv']),'Delimiter',',');
    end
end
if datenum(camerafile(61:end-8),'yymmddHHMMSS') > min(camera.time(4))
    if exist([camerafile(1:end-4) '_ch4_detect.csv'],'file')
        cam_4 = table2array(readtable([camerafile(1:end-4) '_ch4_detect.csv']),'Delimiter',',');
    end
end

%% Load .csv Detection File

detect = table2array(readtable([detectfile(1:end-4) '_check.csv']),'Delimiter',',');

%% Analysis Paramaters

fs = getfield(audioinfo(audiofile),'SampleRate'); % .wav sampling rate (Hz)
array.sample = round(array.tmax.*fs); % maxmimum time delay between channels

%% Check if Detection on Channel 1 is on Channel 2, 3 and 4

disp('Begin check....');

if exist('cam_2','var')
    cam_2_ch1 = [cam_2(:,1)-array.sample(1) cam_2(:,2)+array.sample(1)];
end
if exist('cam_3','var')
    cam_3_ch1 = [cam_3(:,1)-array.sample(2) cam_3(:,2)+array.sample(2)];
end
if exist('cam_4','var')
    cam_4_ch1 = [cam_4(:,1)-array.sample(3) cam_4(:,2)+array.sample(3)];
end

keep = zeros(length(detect),1);
for i = 1:size(detect,1)
    if exist('cam_1','var')
        for j = 1:size(cam_1,1)
            if abs(cam_1(j,1) - detect(i,1)) <= 5*fs
                if abs(cam_1(j,2) - detect(i,2)) <= 5*fs
                    keep(i) = keep(i) + 1;
                end
            end
        end
    end
    if exist('cam_2_ch1','var')
        for j = 1:size(cam_2_ch1,1)
            if abs(cam_2_ch1(j,1) - detect(i,1)) <= 5*fs
                if abs(cam_2_ch1(j,2) - detect(i,2)) <= 5*fs
                    keep(i) = keep(i) + 1;
                end
            end
        end
    end
    if exist('cam_3_ch1','var')
        for j = 1:size(cam_3_ch1,1)
            if abs(cam_3_ch1(j,1) - detect(i,1)) <= 5*fs
                if abs(cam_3_ch1(j,2) - detect(i,2)) <= 5*fs
                    keep(i) = keep(i) + 1;
                end
            end
        end
    end
    if exist('cam_4_ch1','var')
        for j = 1:size(cam_4_ch1,1)
            if abs(cam_4_ch1(j,1) - detect(i,1)) <= 5*fs
                if abs(cam_4_ch1(j,2) - detect(i,2)) <= 5*fs
                    keep(i) = keep(i) + 1;
                end
            end
        end
    end
end

keep  = keep >= 1;

%% Write to Table

if sum(keep) > 0
    var_names = {'StartSample','EndSample'};
    out = table(detect(keep,1),detect(keep,2),'VariableNames',var_names);

    [~,out_file,~]=fileparts(detectfile);
    writetable(out,strcat(out_dir,slash,out_file,'_cam_check.csv'));
end

disp('Finish check!');

end