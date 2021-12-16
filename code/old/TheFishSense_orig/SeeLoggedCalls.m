%% INTRO
% If you have any questions about the code please email me at nsj2@rice.edu

%This can pull calls from the log, the detector, or just randomly.
%For each call this produces 11 figures that walk you through the filtering process that the detector goes through
%This allows you to see how well the filters work with your data

%Sometimes the plots will look like they are getting cut off towards the end, there is no error there
%The plots are the correct length, the axis just continue past the data to the next integer.


clear buffertime timeinterval threshold energypercentile Labels StartSamples EndSamples StartTime EndTime startingnumber timeandindex

%% Getting operating system
%this is so the program is compatible with both pc and mac operating systems
if ispc
    slash='\';
elseif ismac
    slash='/';
else
    disp('Dont recognise operating system');
end

%% What do you want to do
choice = questdlg('What do you want to Look at?', 'Loading Question', 'Manual Log', 'Detections','Random','Detections');

%% Manual Log
if strcmp(choice, 'Manual Log')%for looking at the log
    
    Loadingchoice = questdlg('Select the Log (.xls)', 'Loading', 'OK','Cancel','OK');
    if strcmp(Loadingchoice, 'Cancel')%if you select cancel, program ends
        return
    end
    [LogFile,LogFilepath]=uigetfile('*.xls','Select the Log File (.xls)');%getting the Log
    Log = readtable(fullfile(LogFilepath,LogFile));%importing the log
    Log = Log(:,1:6);%we only care about these columns
    messagebox=msgbox(['These audioFiles have been logged';unique(Log{:,1})],'non-modal');set(messagebox,'Position',[50 400 350 90]);
    
    Loadingchoice = questdlg('Select the .wav or .x.wav file', 'Loading', 'OK','Cancel','OK');
    if strcmp(Loadingchoice, 'Cancel')%if you select cancel, program ends
        return
    end
    
    [audioFile,audioFilepath]=uigetfile('*.wav','Select the .wav or .x.wav File');%getting the audioFile
    delete(messagebox);
    Fs= getfield(audioinfo(fullfile(audioFilepath,audioFile)),'SampleRate');%getting sampling Frequency
    
    
    Log{:,1}=num2cell(transpose(wavname2dnumNic(Log{:,1})));%changing the filename to datenum of the Log
    Log=Log(cell2mat(Log{:,1})==wavname2dnumNic(audioFile),:);%resizing log to only contain the data from the audioFile we are looking at
    
    if isempty(Log)
        disp('No Logged Data for that audioFile, ending program')
        return
    end
    
    StartSamples=round(seconds(duration(Log{:,5}-datetime(cell2mat(Log{1,1}),'ConvertFrom','datenum'),'Format','s')).*Fs);
    EndSamples=round(seconds(duration(Log{:,6}-datetime(cell2mat(Log{1,1}),'ConvertFrom','datenum'),'Format','s')).*Fs);
    Labels=Log{:,4};
    StartTime=Log{:,5};
    EndTime=Log{:,6};
    
    while ~exist('buffertime','var')||isempty(buffertime)||isnan(buffertime)|| sum(buffertime<0);%input scalar value for buffertime. buffertime must be nonnegative
        %while buffertime doesnt exist, while it is empty, while it is NaN, while it is negative
        buffertime=str2double(inputdlg('How many seconds of buffertime would you like on the ends of each call?','Input',1,{'0'}));
    end
    
    while ~exist('threshold','var')||isempty(threshold)||isnan(threshold)|| sum(threshold<0);%input scalar value for threshold. threshold must be nonnegative
        %while threshold doesnt exist, while it is empty, while it is NaN, while it is negative
        threshold=str2double(inputdlg('Input Threshold:','Thresold',1,{'1.5'}));
    end
    
    while ~exist('energypercentile','var')||isempty(energypercentile)||isnan(energypercentile)|| sum(energypercentile<0);%input scalar value for energypercentile. energypercentile must be nonnegative
        %while threshold doesnt exist, while it is empty, while it is NaN, while it is negative
        energypercentile=str2double(inputdlg('Input EnergyPercentile (0-100):','EnergyPercentile',1,{'95'}));
    end
    
end

%% Detections
if strcmp(choice, 'Detections') %for looking at the detections
    
    %Getting audioFile
    Loadingchoice = questdlg('Select the .wav or .x.wav file','Loading', 'OK', 'Cancel','OK');
    if strcmp(Loadingchoice, 'Cancel')%if you select cancel, program ends
        return
    end
    [audioFile,audioFilepath]=uigetfile('*.wav','Select the .wav or .x.wav File');%getting the audioFile
    Fs= getfield(audioinfo(fullfile(audioFilepath,audioFile)),'SampleRate');%getting sampling Frequency
    
    %Getting detectorFile
    Loadingchoice = questdlg('Select the .txt Detector File','Loading', 'OK', 'Cancel','OK');
    if strcmp(Loadingchoice, 'Cancel')%if you select cancel, program ends
        return
    end
    [detectorFile,detectorFilepath]=uigetfile('*.txt','Select the .txt Detector File');%getting .txt file
    underscores=strfind(detectorFile,'_');%finding locations of underscores in detectorFile name
    
    %checking to see if audioFile and detectorFile match
    if wavname2dnumNic(audioFile)~=wavname2dnumNic(detectorFile)%if the audio and detector filenames dont match
        disp(['audioFile = ', audioFile]);
        disp(['detectorFile = ', detectorFile]);
        choosetocontinue = questdlg('audioFile and detectorFile do not appear to correspond, Do you want to continue?', 'Warning!', 'Yes', 'Cancel','Cancel');
        if strcmp(choosetocontinue, 'Cancel')%if you select cancel, program ends
            return
        end
    end
    
    %pulling parameters from the name of the detector file so that these images are processed with the same parameters as the detector
    energypercentile=str2double(detectorFile(underscores(end-1)+1:underscores(end)-1));
    threshold=str2double(detectorFile(underscores(end-3)+1:underscores(end-2)-1));
    
    
    while ~exist('buffertime','var')||isempty(buffertime)||isnan(buffertime)|| sum(buffertime<0);%getting scalar value for buffertime. buffertime must be nonnegative
        %while buffertime doesnt exist, while it is empty, while it is NaN, while it is negative
        buffertime=str2double(inputdlg('How many seconds of buffertime would you like on the ends of each call?','Input',1,{'0'}));
    end
    
    
    testtable=readtable(fullfile(detectorFilepath,detectorFile),'Delimiter',',');%reading the .txt file
    StartSamples=round(table2array(testtable(:,1)).*Fs);%starttime column vector in seconds
    StartSamples(1)=StartSamples(1)+1;%Accounting for the case when StartSamples(1)=0
    EndSamples=round(table2array(testtable(:,2)).*Fs);%endtime column vector in 
    StartTime=datetime(testtable{:,3},'InputFormat','yyyy/MM/dd_HH:mm:ss');%getting starttimes
    EndTime=datetime(testtable{:,4},'InputFormat','yyyy/MM/dd_HH:mm:ss');%getting endtimes
    
    
end



%% Random
if strcmp(choice, 'Random')
    
    %Getting audioFile
    Loadingchoice = questdlg('Select the .wav or .x.wav file', 'Loading','OK', 'Cancel','OK');
    if strcmp(Loadingchoice, 'Cancel')%if you select cancel, program ends
        return
    end
    
    [audioFile,audioFilepath]=uigetfile('*.wav','Select the .wav or .x.wav File');%getting the audioFile
    Fs= getfield(audioinfo(fullfile(audioFilepath,audioFile)),'SampleRate');%getting sampling Frequency
    Duration=getfield(audioinfo(fullfile(audioFilepath,audioFile)),'Duration'); % getting Duration
    
    StartSamples=zeros(200,1);%have to preallocate here so the for loop will run 1:length(StartSamples)
    EndSamples=zeros(200,1);%have to preallocate here
    StartTime=cell(200,1);
    EndTime=cell(200,1);
    
    while ~exist('timeinterval','var')||isempty(timeinterval)||isnan(timeinterval)|| sum(timeinterval<=0);%getting scalar value for timeinterval. must be positive
        %while timeinterval doesnt exist, while it is empty, while it is NaN, while it is not positive
        timeinterval=str2double(inputdlg('How many seconds would you like each random detection to be?','Input',1,{'15'}));
    end
    
    while ~exist('buffertime','var')||isempty(buffertime)||isnan(buffertime)|| sum(buffertime<0); %getting scalar value for buffertime. must be nonnegative4
        %while buffertime doesnt exist, while it is empty, while it is NaN, while it is negative
        buffertime=str2double(inputdlg('How many seconds of buffertime would you like on the ends of each call?','Input',1,{'0'}));
    end
    
    while ~exist('threshold','var')||isempty(threshold)||isnan(threshold)|| sum(threshold<0);%input scalar value for threshold. must be nonnegative
        %while threshold doesnt exist, while it is empty, while it is NaN, while it is negative
        threshold=str2double(inputdlg('Input Threshold:','Thresold',1,{'1.5'}));
    end
    
    while ~exist('energypercentile','var')||isempty(energypercentile)||isnan(energypercentile)|| sum(energypercentile<0);%input scalar value for threshold. must be nonnegative
        %while threshold doesnt exist, while it is empty, while it is NaN, while it is negative
        energypercentile=str2double(inputdlg('Input EnergyPercentile (0-100):','EnergyPercentile',1,{'95'}));
    end
    
end

%% Getting Starting Index
    while ~exist('startingnumber','var')||isempty(startingnumber)||isnan(startingnumber)|| (startingnumber<=0);%input scalar value for startingnumber. startingnumber must be positive
        %while startingnumber doesnt exist, while it is empty, while it is NaN, while it is not positive
        startingnumber=round(str2double(inputdlg(strcat('What call index would you like to start from (1-',string(length(StartTime)),')'),'Input',1,{'1'})));
    end

timeandindex.Position=[554.0000,  475.1667,  152.0000,   74.7500];%preallocating position of timeandindex messagebox

for i= startingnumber:length(StartSamples)
    %% Loading
    
    if strcmp(choice,'Random');%if you clicked 'Random'
        RandomNumber=randi(Duration);
        StartTime{i}=datetime(wavname2dnumNic(audioFile),'ConvertFrom','datenum')+seconds(RandomNumber);%randomly generating starttimes
        EndTime{i}=datetime(wavname2dnumNic(audioFile),'ConvertFrom','datenum')+seconds(RandomNumber+timeinterval);%randomly generating endtimes
        StartSamples(i)=round(RandomNumber*Fs);%getting startsamples
        EndSamples(i)=round(StartSamples(i)+timeinterval*Fs);%getting endsamples
    end
    
    if StartSamples(i)-(buffertime*Fs)<1 %making sure index doesnt excede matrix dimensions due to buffer time
        disp('Skipping Because Index would have exceeded matrix dimensions');
        continue
    end
    
    if EndSamples(i)+(buffertime*Fs)> getfield(audioinfo(fullfile(audioFilepath,audioFile)),'TotalSamples') %making sure index doesnt excede matrix dimensions due to buffer time
        disp('Skipping Because Index would have exceeded matrix dimensions');
        continue
    end
    
    %% Displaying Information about the call
    disp(strcat('StartTime=',string(StartTime(i))));%displays to command window
    disp(strcat('EndTime = ', string(EndTime(i))));%displays to command window
    oldposition=timeandindex.Position;%getting the location of the messagebox( incase people dragged it to a new location)
    timeandindex=msgbox([{strcat('StartTime=',string(StartTime(i)))};{strcat('EndTime = ', string(EndTime(i)))};{strcat('Index = ',string(i))}],'Details','replace');%updating the messagebox
    timeandindex.Position=oldposition;%if you dragged the messagebox to a new location it will update there
    
    %% Pulling the data via audioread
    data=audioread(fullfile(audioFilepath,audioFile),[StartSamples(i)-buffertime*Fs,EndSamples(i)+buffertime*Fs],'native');data=double(data);
    
    %% Original Unfiltered Spectrogram
    figure(1) %Original spectrogram
    [~, frequency, time, psd] = spectrogram(data, hanning(256),158,256,2000,'yaxis','psd');%generating psd, time and frequency
    time=transpose(time);%I like time as a column vector
    imagesc([time(1),time(end)],[frequency(1) frequency(end)],10*log10((psd)));set(gca,'Ydir','Normal');%plotting spectorgram
    c=colorbar; c.Label.String='Power/Frequency (dB/HZ)';%creating colorbar and label
    %caxis([-20 22]);
    colormap parula;
    
    if exist('Labels','var')
        title(strcat(Labels{i},' Original Spectrogram'));% for when we have labels
    else
        title('Original Spectrogram');%for when we dont have labels
    end
    
    xlabel('Time (sec)');
    ylabel('Frequency (Hz)');
    
    
    figure(2) %Original data
    plot(linspace(0,length(data)/Fs,length(data)),data)
    title('Original Data');
    xlabel('Time(seconds)');
    ylabel('Amplitude');
    
    %% Whitening PWelch
    figure(3) %Original PWelch
    plot(frequency,pow2db(sum(psd,2)./length(time)));grid on;
    title('Original Welch Power Spectral Density Estimate');
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    ylim([-10 35]);
    
    figure(4) %Whitened PWelch
    psd=psd./(sum(psd,2)./length(time));%normalizing the pwelch to the average, both methods produce the same results
    %psd=psd./sum(psd,2);%normalizing the pwelch to zero, both methods produce the same results
    plot(frequency,pow2db(sum(psd,2)));grid on;
    title('Whitened Welch Power Spectral Density Estimate');
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    ylim([-10 35]);
    
    %% Whitening Spectrogram and Data
    
    figure(5) %Whitened Spectrogram
    imagesc([time(1),time(end)],[frequency(1) frequency(end)],10*log10((psd)));set(gca,'Ydir','Normal');
    c=colorbar; c.Label.String='Power/Frequency (dB/HZ)';%creating colorbar and label
    % caxis([-15 15]);
    colormap parula
    title('Whitened Spectrogram');
    xlabel('Time (sec)');
    ylabel('Frequency (Hz)');
    
    figure(6) %Whitened Data
    data=transpose(sum(psd));
    plot(time,data);
    title('Whitened Data');
    xlabel('Time (sec)');
    ylim([0 400]);
    ylabel('Power');
    noisetest=mean(findpeaks(data(data>=0)));%finding average noise level
    dispstring={strcat('Mean Noise :',string(noisetest)),strcat('Max Peak to Noise Ratio :',string(max(data)/noisetest))};
    annotation('textbox',[.525 .51 .4 .4],'String',dispstring,'FitBoxToText','on');
    hold on
    plot(time,(threshold*noisetest)*ones(size(data)))
    hold off
    
    %% High-Pass Filter and Boat Filter and DB Filter
    
    figure(7) %Highpass Spectorgram
    psd(1:10,:)=0;%applying the highpass filter
    imagesc([time(1),time(end)],[frequency(1) frequency(end)],10*log10((psd)));set(gca,'Ydir','Normal');
    c=colorbar; c.Label.String='Power/Frequency (dB/HZ)';%creating colorbar and label
    %  caxis([-15 15]);
    colormap parula
    title('Highpass Spectrogram');
    xlabel('Time (sec)');
    ylabel('Frequency (Hz)');
    
    
    
    
    figure(8) %Highpass Data
    data=transpose(sum(psd));
    plot(time,data);
    title('Highpass Data');
    xlabel('Time (sec)');
    ylim([0 400]);
    ylabel('Power');
    noisetest=mean(findpeaks(data(data>=0)));%finding average noise level
    dispstring={strcat('Mean Noise :',string(noisetest)),strcat('Max Peak to Noise Ratio :',string(max(data)/noisetest))};
    annotation('textbox',[.525 .51 .4 .4],'String',dispstring,'FitBoxToText','on');
    hold on
    plot(time,(threshold*noisetest)*ones(size(data)))
    hold off
    
    
    %% DB Filter
    figure(9)%DB Filter Spectrogram
    spectrogramthreshold=prctile(psd(:),energypercentile);%setting the threshold for the db filter
    psd(psd<spectrogramthreshold)=0;%applying the DB Filter
    imagesc([time(1),time(end)],[frequency(1) frequency(end)],10*log10((psd)));set(gca,'Ydir','Normal');
    c=colorbar; c.Label.String='Power/Frequency (dB/HZ)';%creating colorbar and label
    % caxis([5 15]);
    colormap parula
    title('Power Threshold Filter Spectrogram');
    xlabel('Time (sec)');
    ylabel('Frequency (Hz)');
    
    
    figure(10) % DB Filter Data
    data=transpose(sum(psd));
    plot(time,data);
    title('Power Threshold Filter Data');
    xlabel('Time (sec)');
    ylim([0 400]);
    ylabel('Power');
    noisetest=mean(findpeaks(data(data>=0)));%finding average noise level
    dispstring={strcat('Mean Noise :',string(noisetest)),strcat('Max Peak to Noise Ratio :',string(max(data)/noisetest))};
    annotation('textbox',[.525 .51 .4 .4],'String',dispstring,'FitBoxToText','on');
    hold on
    plot(time,(threshold*noisetest)*ones(size(data)))
    hold off
    
    %% Find Peaks Filter
    figure(11) %
    findpeaks(data,time,'MinPeakHeight',threshold*noisetest);%finding the peaks that are larger than the noise
    title('Time Series Threshold Filter Data')
    xlabel('Time (sec)');
    ylabel('Power');
    ylim([0 400]);
    dispstring={strcat('Mean Noise :',string(noisetest)),strcat('Max Peak to Noise Ratio :',string(max(data)/noisetest))};
    annotation('textbox',[.525 .51 .4 .4],'String',dispstring,'FitBoxToText','on');
    hold on
    plot(time,(threshold*noisetest)*ones(size(data)));
    hold off
    
    
    k=waitforbuttonpress;
    delete(findall(figure(6),'type','annotation'));
    delete(findall(figure(8),'type','annotation'));
    delete(findall(figure(10),'type','annotation'));
    delete(findall(figure(11),'type','annotation'));
    
end



