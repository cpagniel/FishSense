function TheFishSense(app)

%% Check Operating System
if ispc % PC
    slash='\';
elseif ismac % Mac
    slash='/';
else
    error('Invalid OS');
end

%% Obtain Audio Files (.wav) of the Input Directory
files = dir(strcat(app.InputDirectoryLocationLabel.Text,slash,'*.wav'));

%% Total Duration of Data in Input Directory
d = zeros(size(files,1),1); % pre-allocating
for i = 1:length(files)
    d(i) = getfield(audioinfo(fullfile(files(i).folder,files(i).name)),'Duration');
end

%% Create Progress Bar for Detection
progress_bar = waitbar(0,['Total Duration of Deployment = ',...
    num2str(round(sum(d)/(60*60*24))),' days']);
pause(1);
start_t = tic;

waitbar(0,progress_bar,'Initializing TheFishSense Detector...');
for i = 1:length(files)
    
    %% Pre-Allocation
    T = [];
    
    %% Audio File Information
    audiofile = fullfile(files(i).folder,files(i).name);
    fs = getfield(audioinfo(audiofile),'SampleRate'); % .wav sampling rate (Hz)
    
    %% DAQ Information
    array = load(app.DAQFileLabel.Text);
    hsens = array(:,1); % hydrophone sensitivity
    tf = array(:,2); % transfer function for DAQ
    
    %% Sample
    
    start_sample = transpose(1:app.SampleDuration.Value*fs:getfield(audioinfo(audiofile),...
        'TotalSamples')); % start sample of each sub sample
    end_sample = min(start_sample + app.SampleDuration.Value*fs - 1, getfield(audioinfo(audiofile),...
        'TotalSamples')); % end sample of each sub sample
    
    for j = 1:(length(start_sample)-1)
        
        %% Detector
        if i == 1
            waitbar(0,progress_bar,'Starting TheFishSense Detector...');
        end
        t = detector(app,audiofile,fs,start_sample(j),end_sample(j),hsens,tf);
        
        %% Correct Times for Sub Sample Start
        for k = 1:size(t,1)
            t(k,:) = (start_sample(j)/fs)+t(k,:);
        end
        
        %% Concatenate
        T = vertcat(T,t);
        clear t
        
    end
    
    %% Write to Table
    
    if ~isempty(T)
        var_names = {'StartSample','EndSample'};
        out = table(T(:,1)*fs,T(:,2)*fs,'VariableNames',var_names);
        [~,out_file,~] = fileparts(audiofile);
        writetable(out,strcat(app.OutputDirectoryLocationLabel.Text,slash,out_file,'_detect.csv'));
    end
    
    %% Update Progress Bar
    elapsed_t = toc(start_t)/(60);
    ETA = ((elapsed_t/sum(d(1:i)))*(sum(d)-sum(d(1:i))));
    waitbar((sum(d(1:i))/sum(d)),progress_bar,sprintf('Overall Progress: %d%%\n ETA: %d min',...
        round((sum(d(1:i))/sum(d))*100),round(ETA)));
end

end