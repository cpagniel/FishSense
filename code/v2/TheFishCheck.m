function TheFishCheck(audiofile,detectfile,out_dir,array)
%% The Fish Check
% Check if detection is on all four channels
%
% INPUTS:
%   audiofile - string of complete path and name of .wav file (e.g., F:Detector\otherstuff\name.wav)
%   detectfile - string of complete path and name of .csv file with initial
%   detections
%   out_dir - string of output directory to save .txt file
%   array - 
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

%% Load .csv Detection File

detect_1 = table2array(readtable([detectfile(1:end-4) '_ch1_detect.csv']),'Delimiter',',');
detect_2 = table2array(readtable([detectfile(1:end-4) '_ch2_detect.csv']),'Delimiter',',');
detect_3 = table2array(readtable([detectfile(1:end-4) '_ch3_detect.csv']),'Delimiter',',');
detect_4 = table2array(readtable([detectfile(1:end-4) '_ch4_detect.csv']),'Delimiter',',');

%% Analysis Paramaters

fs = getfield(audioinfo(audiofile),'SampleRate'); % .wav sampling rate (Hz)
array.sample = round(array.tmax.*fs); % maxmimum time delay between channels

%% Check if Detection on Channel 1 is on Channel 2, 3 and 4

disp('Begin check....');

detect_1_ch2 = [detect_1(:,1)-array.sample(1) detect_1(:,2)+array.sample(1)];
detect_1_ch3 = [detect_1(:,1)-array.sample(2) detect_1(:,2)+array.sample(2)];
detect_1_ch4 = [detect_1(:,1)-array.sample(3) detect_1(:,2)+array.sample(3)];

keep = zeros(length(detect_1),1);
for i = 1:size(detect_1,1)
    for j = 1:size(detect_2,1)
        if detect_2(j,1) >= detect_1_ch2(i,1)
            if detect_2(j,2) <= detect_1_ch2(i,2)
                keep(i) = keep(i) + 1;
            end
        end
    end
    for j = 1:size(detect_3,1)
        if detect_3(j,1) >= detect_1_ch3(i,1)
            if detect_3(j,2) <= detect_1_ch3(i,2)
                keep(i) = keep(i) + 1;
            end
        end
    end 
    for j = 1:size(detect_4,1)
        if detect_4(j,1) >= detect_1_ch4(i,1)
            if detect_4(j,2) <= detect_1_ch4(i,2)
                keep(i) = keep(i) + 1;
            end
        end
    end 
end

keep  = keep == 4;

%% Write to Table

var_names = {'StartSample','EndSample'};
out = table(detect_1(keep,1),detect_1(keep,2),'VariableNames',var_names);

[~,out_file,~]=fileparts(detectfile);
writetable(out,strcat(out_dir,slash,out_file,'_check.csv'));

disp('Finish check!');

end