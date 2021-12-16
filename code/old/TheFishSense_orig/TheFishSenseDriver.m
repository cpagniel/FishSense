%% INTRO
% If you have any questions about the code please email me at nsj2@rice.edu
%  WORKS WITH MULTIPLE CHANNELS- makes a new txt file for each channel
%This CruchTimeDriver can process audioFiles from multiple folders across multiple locations and save the outputted txt files in one location or in
%the same location as the audioFile

%This works by first prompting the user to get the directories of the audioFiles that need to be processed, and also prompts the user for the location
%to save the .txt files

%after that it just runs TheFishSense on all of the files, and tells TheFishSense where to save the .txt files
%TheFishSense is a 2 stage threshold detector

%% Getting operating system
%this is so the program is compatible with both pc and mac operating systems
if ispc
    slash='\';
elseif ismac
    slash='/';
else
    disp('Dont recognise operating system');
end

%% Getting first audioFiles Directory
choice = questdlg('Select the directory of the audiofiles', 'Select Directory', 'OK', 'Cancel','OK');

if strcmp(choice, 'Cancel')%if cancel program ends
    return
end
audiofilesdir=uigetdir;%The folder where the audiofiles are

a = dir(strcat(audiofilesdir,slash,'*.wav'));%structure that contains information about the files to be processed
disp('Folders that have been added');%shows the folders that have been added to the audioFiles Directory
disp(unique({a(:).folder}'));%shows the folders that have been added to the audioFiles Directory
%% Getting additional audioFiles Directory
choice='YES';
while strcmp(choice,'YES')%while choice equals yes
    choice = questdlg('Do you want to add additional directories to process?', 'Select Another Directory', 'YES', 'NO','YES');
    
    if strcmp(choice,'YES')
        audiofilesdir=uigetdir;%The place where the audiofiles are
        a = [a;dir(strcat(audiofilesdir,slash,'*.wav'))];%a grows with each addiditonal directory that is added
        disp('Folders that have been added');%shows the folders that have been added to the audioFiles Directory
        disp(unique({a(:).folder}'));%shows the folders that have been added to the audioFiles Directory
    end
end

%% Getting Save Directory

savedirchoice = questdlg('Do you want to save the .txt files all in one location or do you want the .txt file to be in the same folder as the audiofile?', 'Save Directory', 'One Location', 'With audioFile','One Location');

if strcmp(savedirchoice,'One Location')%if you answered one location then select the location you want to save the .txt files to
    choice = questdlg('Select a directory to store .txt files', 'Select Directory', 'OK', 'Cancel','OK');
    
    if strcmp(choice, 'Cancel')%if cancel program ends
        return
    end
    
    savedir=uigetdir;%place where you want the .txt files to be saved
end

%% Duration Calculations for super cool eta/waitbar
Duration=zeros(size(a,1),1);%preallocating
for i=1:length(a)%iterating through the files that were in the folders of the audioFile directory
    Duration(i)=getfield(audioinfo(fullfile(a(i).folder,a(i).name)),'Duration');%Duration in seconds
end
Duration=Duration./(60);%Duration in min
disp(['Total Data Duration =',' ',num2str(round(sum(Duration)/(60*24))),' ','days']);%Total duration of all the files to be processed

%% Preallocating
errors=[];

%% Doing the Stuff
progressBar = waitbar(0, 'Processing Data Set...','Position',[500 500 275 60]);
tstart1=tic;%important to use tstart1 so the tic here and the tic in TheFishSense dont interact
for i=1:length(a)%iterating through the files that were in the folders of the audioFile directory
    %[intervals] = TheFishSense(audioFile, savedir, channel, chunkDuration, threshold, intervalLength, energypercentile)
    
    for j=1:getfield(audioinfo(fullfile(a(i).folder,a(i).name)),'NumChannels');%iterating through the number of channels

        if strcmp(savedirchoice,'One Location')%savedir type
            [~,errortimes]=TheFishSense(fullfile(a(i).folder,a(i).name), savedir, j , 15, 1.5, 1.5, 95);
        else %if you want the.txt file to save in the same location as the audioFile
            [~,errortimes]=TheFishSense(fullfile(a(i).folder,a(i).name), a(i).folder, j , 15, 1.5, 1.5, 95);
        end

    end

    %Dealing with errors in recordings
    if ~isempty(errortimes)
        errortimes=cellstr(errortimes);
        errortimes(:,2)={fullfile(a(i).folder,a(i).name)};
        errors=vertcat(errors,errortimes);
    end
    
    %waitbar
    elapsedtime=toc(tstart1)/(60);%elapsedtime in min
    eta=((elapsedtime/sum(Duration(1:i)))*(sum(Duration)-sum(Duration(1:i))));%Eta in whatever units the elapsed time and duration are in
    waitbar((sum(Duration(1:i))/sum(Duration)), progressBar, sprintf('Processing Data Set Progress: %d%%\n ETA: %d min',round((sum(Duration(1:i))/sum(Duration))*100),round(eta)));%updating progressbar
end
close(progressBar);
disp('Done!');
