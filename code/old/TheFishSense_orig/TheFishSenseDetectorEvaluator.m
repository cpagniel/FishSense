%% INTRO
% If you have any questions about the code please email me at nsj2@rice.edu

%Used to evaluate the detector
%saves cell array DetectorResults in the current folder
%DetectorResults contains all of the results from the calculator
%Need a manual log and a detector file

%INPUTS:
%   detectorfilesdir - user is prompted to input
%   Log - user is prompted to input
%
%OUTPUTS:
%   ChunkDuration - detector parameter
%   Threshold - detector parameter
%   Interval - detector parameter
%   Energy_Percentile - detector parameter
%   Total_Points - total number of detections that detector produced
%   Duration - median duration of detections 
%   Accuracy - percentage of logged calls that the detector was able to detect
%   MultipleCall% - measures a percentage of how many of the detections had more than one logged call in that single detection, see picture
%   Intersections- column vector containing times of each logged call that the detector picked up
%   Number - histogram breakdown of how many logged calls were present in each detection
%   Type - histogram breakdown of the type of calls that the detector picks up
%   Filename - the name of each detectorfile just to be sure everything lines up

warning off
progressBar = waitbar(0, 'Processing variations...');%making progressBar

%% Getting operating system
%this is so the program is compatible with both pc and mac operating systems
if ispc
    slash='\';
elseif ismac
    slash='/';
else
    disp('Dont recognise operating system');
end

%% Getting Detector .txt files directory
choice = questdlg('Select the directory of the detector .txt files', 'Select Directory', 'OK', 'Cancel','OK');

if strcmp(choice, 'Cancel')
    return
end
detectorfilesdir=uigetdir;%The folder where the audiofiles are
a = dir(strcat(detectorfilesdir,slash,'*.txt'));%structure that contains information about the detectorfiles to be processed
%% Getting the Log
Loadingchoice = questdlg('Select the Manual Log (.xls)', 'Loading', 'OK','Cancel','OK');
if strcmp(Loadingchoice, 'Cancel')%if you select cancel, program ends
    return
end
[LogFile,LogFilepath]=uigetfile('*.xls','Select the Log File (.xls)');%getting the Log
Log=readtable(fullfile(LogFilepath,LogFile));%importing the log
Log=Log(:,1:6);%only care about these columns in the log

Log{:,1}=num2cell(transpose(wavname2dnumNic(Log{:,1})));%changing the filename to datenum in the Log
a=a(ismember(wavname2dnumNic({a.name}),unique(cell2mat(table2array(Log(:,1))))));%removing the files that havent been logged
%% Preallocating
DetectorResults=cell(size(a,1),12); %preallocating cell size
%I keep the Column Headers in a format to also be variable names so I can easily convert to a table when I want to (thats why I use the underscores)
DetectorResults{1,1}='Chunk_Duration';DetectorResults{1,2}='Threshold';DetectorResults{1,3}='Interval';DetectorResults{1,4}='Energy_Percentile';DetectorResults{1,5}='Total_points';DetectorResults{1,6}='Duration';DetectorResults{1,7}='Accuracy';DetectorResults{1,8}='MultipleCall%';DetectorResults{1,9}='Intersections';DetectorResults{1,10}='Number';DetectorResults{1,11}='Type';DetectorResults{1,12}='Filename';


tic; %start timer for eta calculations
for i=1:size(a,1)
    %% Adjusting the log for the audioFile of the Detector
    LogNew=Log(cell2mat(Log{:,1})==wavname2dnumNic(fullfile(a(i).folder,a(i).name)),:);%resizing log to only contain the data from the audioFile we are looking at
    
    if isempty(LogNew)%if there is no logged data for that detector file
        disp('No Logged Data for that detector file, skipping it');
        disp(a(i).name);
        continue
    end
    Labels=LogNew{:,4};
    startTimereal=LogNew{:,5};
    endTimereal=LogNew{:,6};
    
    %Removing repeats from the log
    totaltest=strcat(datestr(startTimereal),datestr(endTimereal),Labels);%only removing exact repetitions in starttime,endtime,and labels
    %cant just use unique(startTimereal),unique(endTimereal),unique(Labels), because after the indicies will no longer correspond to each other
    [~,e,~]=intersect(totaltest,unique(totaltest),'rows');%e is a logical, that specifies the unique rows
    Labels=Labels(e);
    startTimereal=datenum(startTimereal(e));
    endTimereal=datenum(endTimereal(e));
    
    
    %% Calculations
    
    filename=a(i).name;%loading one txt file
    underscores=strfind(strrep(filename,'Filter.txt',''),'_');%pulling the parameters from its name
    variables=filename((underscores(end-4))+1:(underscores(end))-1);%pulling the parameters from its name
    
    testtable=readtable(fullfile(a(i).folder,a(i).name),'Delimiter',',');%reading the txt file
    startduration=table2array(testtable(:,1));%starttime in seconds with respect to the start of the recording
    endduration=table2array(testtable(:,2));%endtime in seconds with respect to the start of the recording
    startTimetest=datenum(table2array(testtable(:,3)));%global starttime
    endTimetest=datenum(table2array(testtable(:,4)));%global endtime
    
%     startTimetest=startTimetest((endduration-startduration)<30);%Getting rid of intervals that are ridiculously long
%     endTimetest=endTimetest((endduration-startduration)<30);%Getting rid of intervals that are ridiculously long
    
    intersections=cell(length(startTimereal),1);%preallocating
    number=zeros(length(startTimereal),1);%preallocating
    
    %checking if any of the logged times fit within the detection times
    for j=1:length(startTimetest)%iterating through all of the detection windows
        tlower= startTimetest(j);
        tupper= endTimetest(j);
        %I am adding a 2 second buffer to the detection window to account for discripancies in the log
        intersections{j}= datestr(startTimereal(startTimereal >= addtodate(tlower,-2,'second') & endTimereal <= addtodate(tupper,2,'second')),'yyyy/mm/dd_HH:MM:SS.FFF');%all of the intersections
        number(j)=size(intersections{j},1);%how many logged calls fit in that detection window
    end
    intersections = intersections(~cellfun('isempty',intersections));%removing the emptys from intersections
    
    %Making intersections one long list of datestrs,
    b=[];%preallocating
    for j=1:length(intersections)%sometimes an intersection cell has more than one datestr, so we have to use this for loop to iterate through the contents of each intersections cell
        b=[b;intersections{j}];
    end
    intersections=cellstr(b);%now intersections is one long list of datestrings
    
    [counts,bins]=histcounts(categorical(number));%getting histogram breakdown of number
    underscores=strfind(variables,'_');%pulling parameters from name of file
    
    DetectorResults{i+1,1}= str2double(variables(1:(underscores(1)-1)));%Chunk Duration(seconds)
    DetectorResults{i+1,2}= str2double(variables(underscores(1)+1:underscores(2)-1));%Threshold
    DetectorResults{i+1,3}= str2double(variables(underscores(2)+1:underscores(3)-1));%Interval
    DetectorResults{i+1,4}= str2double(variables(underscores(3)+1:end));%Energy Percentile
    DetectorResults{i+1,5}= length(startTimetest);%Total_Points
    DetectorResults{i+1,6}= median(endduration-startduration);%Median Duration
    DetectorResults{i+1,7}=(sum(number)/length(startTimereal))*100;%Accuracy, percent of logged calls that were detected(does not remove multiple calls)
    DetectorResults{i+1,8}=(sum(counts(3:end))/(sum(counts(2:end))))*100;%Multiple Call Percentage, 
    DetectorResults{i+1,9}= intersections;%Intersections
    DetectorResults{i+1,10}=[transpose(bins),num2cell(transpose(counts))];%Histogram break down of number
    
    Index=[];%finding the indicies of logged calls that the detector was able to detect
    for j=1:length(intersections)
        Index=[Index; find(contains(datestr(startTimereal,'yyyy/mm/dd_HH:MM:SS.FFF'),intersections(j,:)))];
    end
 %% Getting the type of each logged call that was detected
 %Getting the tables to append to eachother based on the Type was a pain
 
    %Breakdown of Labels from the log
    [Typecounts,Type]=histcounts(categorical(Labels));%Type is the label, Typecounts are the times that label appeared
    Type=transpose(Type);%making column vector
    Typecounts=transpose(Typecounts);%making column vector
    TypeBreakdown=table(Type,Typecounts);%creating table
    TypeBreakdown.Properties.VariableNames{2}='Log';%first column of table is the variable name Type, second column of tabel is variable name Log
    
    %Breakdown with multiple calls
    [Typecounts,Type]=histcounts(categorical(Labels(Index)));%looking at all intersections
    Typecounts=transpose(Typecounts);Type=transpose(Type);%making column vector
    TypeBreakdownNew=table(Type,Typecounts);%making a new table
    TypeBreakdownNew.Properties.VariableNames{2}='Repeats';%first column of table is the variable name Type, second column of tabel is variable name Repeats
    TypeBreakdown=outerjoin(TypeBreakdown,TypeBreakdownNew,'MergeKeys',true);%the key variable for the merge is Type, 
    %now TypeBreakdown is a 3 column table: Type,Log,Repeats
    
    %Breakdown with only unique calls
    [Typecounts,Type]=histcounts(categorical(Labels(unique(Index))));%only looking at the unique indicies of intersections
    Typecounts=transpose(Typecounts);Type=transpose(Type);%making column vector
    TypeBreakdownNew=table(Type,Typecounts);%making new table
    TypeBreakdownNew.Properties.VariableNames{2}='Unique';%first column of table is the variable name Type, second column of tabel is variable name Unique
    TypeBreakdown=outerjoin(TypeBreakdown,TypeBreakdownNew,'MergeKeys',true);%key variable for the merge is Type
    %now TypeBreakdown is a 4 column table: Type, Log, Repeats, Unique
    
    
    DetectorResults{i+1,11}=TypeBreakdown;%Histogram of Types
    DetectorResults{i+1,12}=filename;%filename
    
    elapsedtime=toc;
    eta=((elapsedtime/(i))*(size(a,1)-(i)))/60;%calculating eta
    waitbar(i/(size(a,1)), progressBar, sprintf('Processing Progress: %d%%\n ETA: %d min', round(i/(size(a,1))*100),round(eta)));%updating progressBar
    
end

save ('DetectorWeekend', 'DetectorResults');%saving in current folder
%% To Make and Write to Excel Table
%DetectorResultstable=array2table(horzcat(DetectorResults(:,1:8),DetectorResults(:,12)));
%writetable(DetectorResultstable,'F:\Detector\FinalReport\DetectorFinalReportResults.xlsx','Sheet','Sheet1');%you might need to change save location

%% Done
close(progressBar);
disp('Done!');

