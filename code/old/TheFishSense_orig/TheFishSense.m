function [intervals,errortimes] = TheFishSense(audioFile, savedir, channel, chunkDuration, threshold, intervalLength, energypercentile)
%% Intro
% This detects more than fish, it detects loud sounds. I was going to call it the LoudSoundDetector but idk about those abreviations.
% If you want to look at previous versions of this code, it used to be named mainFiltered
% If you have any questions about the code please email me at nsj2@rice.edu
%
%
%INPUTS:
%   audioFile- string that is the complete path and name of the file (EX: F:Detector\otherstuff\name.x.wav (or just.wav))
%   savedir - string that is the location of where you want to save the .txt file
%   channel - specifies which channel to look at
%   chunkDuration- How long each chunk is (in seconds)
%   threshold- a scalar that you multiply by the mean noise level to get the threshold level
%   intervalLength- scalar in units of seconds, used to make start and end times from the peak locations
%                   Also key because it affects how long there must be no peaks between calls
%   energypercentile- a scalar (0-100) that is used to set the threshold of the power filter
%
% WARNING DON'T CHOOSE A CHUNKDURATION TOO BIG OR YOU WILL OVERLOAD THE RAM ON THE COMPUTER AND IT WILL CRASH
% SPEAKING FROM EXPERIENCE
% TRYING TO READ THE ENTIRE AUDIOFILE WILL CRASH THE COMPUTER
%
%OUTPUTS:
%   Intervals: a table that contains 4 columns, the first column in starttimes in seconds, second: endtimes in seconds,third:starttimes global, fourth:endtimes global
%   Intervals is written into a .txt document and saved to the location specified by savedir
%   errortimes- the times that there appears to be an error in the hydrophone recording

%% Getting operating system
%this is so the program is compatible with both pc and mac operating systems
if ispc
    slash='\';
elseif ismac
    slash='/';
else
    disp('Dont recognise operating system');
end
%% Loading stuff

tic %starting timer for eta measurements
disp('Program starts....');
result = [];%preallocating
errortimes={};%preallocating
globaltime=wavname2dnumNic(audioFile);%using this function to pull the global starttime of the file, in datenum format

progressBar = waitbar(0, 'Processing audiofile...');%making progressBar
Fs = getfield(audioinfo(audioFile),'SampleRate');%getting sampling rate of audiofile

%% Chunking the audio segments and analyzing the chunks
disp('Start to analyse signal and find peaks....');


Chunks = transpose([1:chunkDuration*Fs:getfield(audioinfo(audioFile),'TotalSamples')]);%building a vector that contains the start chunk samples
for i=1:(length(Chunks)-1)%iterating through each start sample in the Chunks vector
    startLoc = Chunks(i);%start location is the sample in chunks
    endLoc = min(startLoc + chunkDuration*Fs - 1, getfield(audioinfo(audioFile),'TotalSamples'));%end location is either the sample before the next start location or the end of the file
    
    [samples,~] = audioread(audioFile, [startLoc endLoc],'native');%reading the audiofile native works best
    samples=samples(:,channel);%only looking at the data of the channel we are looking at
    if Chunks(i)==Chunks(end)
        2+2;
    end

    elapsedtime=toc/60;%elapsedtime in minutes
    eta=((elapsedtime/(i))*(size(Chunks,1)-(i)));%eta in whatever units elapsedtime is in
    waitbar(i/(size(Chunks,1)), progressBar, sprintf('Processing audiofile progress: %d%%\n ETA: %d min', round(i/(size(Chunks,1))*100),round(eta)));%upadating progressBar
    
    [signals, errorflag] = signal_detectionFiltered(samples, Fs, chunkDuration, threshold, energypercentile);%Where the work gets done
    signals = (signals + (startLoc-1)/Fs); %shifting the times
    result = vertcat(result, signals);%appending the times
    
    if errorflag==1 %if there was an error in this chunk, errorflag will be 1
        errortimes=vertcat(errortimes,datetime(addtodate(globaltime,Chunks(i)/Fs,'second'),'ConvertFrom','datenum'));%the chunk that contains the error
    end
    
end

disp('Finish!');
%% Creating Intervals Without Overlap

%This is where we turn the peak times into start and end times
%This method produces the same results as the chaining method, but much faster and cleaner
%Because this method does matrix operations instead of going through the detections with scalar operations
disp('Start to create intervals and merge them....');

difference=find(diff(result)>2*intervalLength);%finding the locations that have more than 2*interval length inbetween them
startTimesort=[result(1)-intervalLength;result(difference(:)+1)-intervalLength];%making start times
startTimesort(startTimesort<0)=0;%if the starttimesort is less than zero, make it zero

endTimesort=[result(difference(:))+intervalLength;result(end)+intervalLength];% making end times
endTimesort(endTimesort>getfield(audioinfo(audioFile),'Duration'))= getfield(audioinfo(audioFile),'Duration');%if endtime is after the end of the file, make it the end of the file


%Now startTimesort and endTimesort are our start and end times without overlaps
disp('Finish!');

%% Making the starttimes and endtimes global

%Until now the start and end times have been seconds in reference to the start of the file
%This converts those times into date strings that give year/month/day_Hour(24)/minute/second
disp('Making times global');
startTimeg=cell(length(startTimesort),1);endTimeg=cell(length(startTimesort),1);%preallocating

%We are taking our global starttime of the file, and adding our startTimesorts to that so we get global starttimes
%Then we do the samething with endtimes
for i=1:length(startTimesort)%iterating through our startTimesort
    startTimeg{i,1}=datestr(addtodate(globaltime,round(startTimesort(i)),'second'),'yyyy/mm/dd_HH:MM:SS');%adding each startTimesort to the globalstart time
    endTimeg{i,1}=datestr(addtodate(globaltime,round(endTimesort(i)),'second'),'yyyy/mm/dd_HH:MM:SS');%adding each endTimesort to the global endtime
end
disp('Finish!');

%now we have 4 column vectors that are all the same size
%startTimesort =starttimes in seconds (column vector matrix)
%endTimesor t= endtimes in seconds (column vector matrix)
%startTimeg = starttimes in global time (column vector cell)
%endTimeg = endtimes in global time (column vector cell)
%% Exporting to txt
disp('Exporting to .txt');

[~,audioFilename,~]=fileparts(audioFile);%Since our audioFile is the entire path to the file. We need to break that down to get only the name of the file.

outputFileName = sprintf('%s_%s_intervals_%d_%g_%g_%g_', audioFilename,strcat('CH',string(channel)), chunkDuration, threshold, intervalLength, energypercentile);%building output filename based on the name of the file and parameters used
intervals=array2table([num2cell(startTimesort),num2cell(endTimesort),startTimeg,endTimeg]);%making table (4 columns startTimesort, endTimesort, startTimeg, endTimeg)

writetable(intervals,strcat(savedir,slash,outputFileName,'.txt'),'WriteVariableNames',0);%writing the table to the location specified by savdir

if ~isempty(errortimes)%if there was an error, make a table containing the errors
    writetable(table(errortimes),strcat(savedir,slash,outputFileName,'Errors.txt'),'WriteVariableNames',1);%wrtiting table of errors to the location
end

disp('Done!');
close(progressBar);
toc
end

function [ locations,errorflag ] = signal_detectionFiltered( data, frequency, ~, threshold, energypercentile)
%% INTRO
% If you have any questions about the code please email me at nsj2@rice.edu
%INPUTS:
%   data- a vector containing the samples from an audio file, comes in as int16 format
%   frequency- the sampling frequency from the audiofile
%   chunkDuration- the scalar time in seconds that the file should be restricted to (~ because we dont need it anymore)
%   threshold- a scalar that you multiply by the mean noise level to get the threshold level
%   energypercentile- a scalar (0-100) that is used to set the threshold of the power filter

%OUTPUTS: column vector of peakTimes for all the hits in this audiofile
%% Preallocating
errorflag=0;%
%% Filtering
%I know this method of Filtering is unorthodox, but it works better for what we are doing than using a highpass filter in time series

data=double(data);%very important, needs to be double format

%Getting psd to do filtering with the matrix
[~, ~, time, psd] = spectrogram(data, hanning(256),128,256,frequency,'yaxis','psd');
time=transpose(time);%making time a column vector

%Whitening by Normalizing Noise with respect to frequency
%if the sum(psd,2) is zero then the psd becomes NaN when dividing by zero
psd=psd./sum(psd,2);%normalizing to zero, both methods get the same result
%psd=psd./(sum(psd,2)./length(time));%normalizing the pwelch to the average noise level, both methods get the same result

%HighPass Filter
psd(1:10,:)=0;%applying the highpass filter

%First Power Threshold
spectrogramthreshold=prctile(psd(:),energypercentile);%setting the threshold for the db filter
psd(psd<spectrogramthreshold)=0;%applying the DB Filter

%Summing with respect to time, getting back to spectra
if ~isempty(find(isnan(psd),1)) %if there are nan values in the psd matrix
    data=transpose(nansum(psd));%Using nansum bc sometimes the amplitude is weird and all zeros, which means my psd is all NaN's
    errorflag=1;
else
    data=transpose(sum(psd));%Summing with respect to time, getting back to spectra
end

%data values need to be in a column vector because of how they are appended in mainFiltered.m
%% Final Threshold

% Find peaks dataLength
noise=mean(findpeaks(data(data>=0)));%finding average noise level
[~, locations] = findpeaks(data,time,'MinPeakHeight',threshold*noise);%locations of peaks are in units of time(in seconds from the start of the audiofile)
end
