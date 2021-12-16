%% The Fish Sense Review

%% Check Operating System

if ispc % PC
    slash='\';
elseif ismac % Mac
    slash='/';
else
    error('Invalid OS');
end

%% Hydrophone Calibration

array.file = 'F:\C_files\C_Winter_2020\july2018_hydrophone_locations.txt'; % file containing array geometry and hydrophone sensitivity
array_raw = load(array.file);

array.hsens = array_raw(:,4); % hydrophone sensitivity for channels 1 to 4 (from calibration sheet)
array.tf = 20*log10(2^16/(1.5--1.5)); % transfer function for ST4300 DAQ

clear array_raw

%% Wav File Diretory

wav_dir = uigetdir(pwd,'Select directory of .wav files');

%% Detection File Directory

f_dir = uigetdir(pwd,'Select directory of .csv files');
files = dir(strcat(f_dir,slash,'*.csv')); % structure that contains information about the files to be processed

%% Additional Directories

choice='YES';
while strcmp(choice,'YES')
    choice = questdlg('Do you want to add additional directories to process?','Select Another Directory','YES','NO','YES');
    if strcmp(choice,'YES')
        f_dir = uigetdir(pwd,'Select additional directory of .csv files');
        files = [files; dir(strcat(f_dir,slash,'*.csv'))];
    end
end

%% Output Directory

out_dir = uigetdir(pwd,'Select directory to save output .csv files');

numb_detect = zeros(size(files,1),1);
numb_detect_true = zeros(size(files,1),1);

%% Review

for i = 6:length(files)
    
    %% Select and load output .csv file location
    
    detect = table2array(readtable(fullfile(files(i).folder,files(i).name),'Delimiter',','));
    if ~isempty(detect)
        numb_detect(i) = size(detect,1);
        
        detect_final = detect;
        
        %% Convert Date/Time to Samples
        
        audiofile = [fullfile(wav_dir,files(i).name(1:end-15)) '.wav'];
        fs = getfield(audioinfo(audiofile),'SampleRate'); % .wav sampling rate (Hz)
        
        % Global Time
        global_t = wavname2dnum_edit(audiofile); % using this function to pull the global starttime of the file, in datenum format
        
        %% Extraction, Plotting and Review
        
        extract = [floor(detect(:,1))-fs ceil(detect(:,2))+fs]; % extract +/- 1 second before and after call
        extract(extract < 0) = 1;
        
        %% Analysis Parameters
        
        channel = [1 2 3 4];
        Nfft = 512; % FFT length
        overlap = 0.9; % 90% overlap
        
        for j = 1:(size(extract,1))
            
            fig = figure();
            movegui('west');
            
            for k = 1:length(channel)
                [data,~] = audioread(audiofile, [extract(j,1) extract(j,2)],'native'); % load .wav file sub sample
                data = data(:,channel(k));
                data = double(data);
                data = 10.^(-1*(array.hsens(k)+array.tf)/20)*data; % DAQ calibration
                data = data-mean(data); % remove mean
                
                %% Spectrogram
                
                subplot(2,2,channel(k))
                [~,F,T,P] = spectrogram(data,kaiser(Nfft,7.85),round(Nfft*overlap),Nfft,fs,'yaxis');
                T = transpose(T);
                
                imagesc(T,F,10*log10(P))
                set(gca,'YDir','normal'); xlim([0 max(T)]); ylim([0 1500]);
                xlabel('Time (s)'); ylabel('Frequency (Hz)'); title([num2str(j) ' of ' num2str(length(extract)) ' detections']);
                colormap(jet); h = colorbar; clim = caxis; caxis([60 90]); % caxis([round(clim(2))-40 clim(2)]);
                ylabel(h,{'Spectral Density (dB re counts^2/Hz)'});
                
                hold on
                plot(T(1)+1.*ones(2,1),[0 1500],'w-','linewidth',2)
                plot(T(end)-1.*ones(2,1),[0 1500],'w-','linewidth',2)
                
            end
            
            %% Review
            
            answer = questdlg('Is this a fish call?',...
                'Keep or Discard Detection','Yes','No','No');
            if strcmp(answer,'No')
                detect_final(j,:) = NaN;
            end
            
            close(fig)
            
        end
        
        detect_final(isnan(detect_final(:,1)),:) = [];
        numb_detect_true(i) = size(detect_final,1);
        
        %% Write to Table
        
        var_names = {'StartSample','EndSample','StartTime','EndTime'};
        
        start_time = cell(size(detect_final,1),1);
        end_time = cell(size(detect_final,1),1);
        for j = 1:size(detect_final,1)
            start_time{j} = datestr(addtodate(global_t,round(detect_final(j,1)./fs*1000),'millisecond'),'yyyy/mm/dd HH:MM:SS:FFF');
            end_time{j} = datestr(addtodate(global_t,round(detect_final(j,2)./fs*1000),'millisecond'),'yyyy/mm/dd HH:MM:SS:FFF');
        end
        
        out = table(detect_final(:,1),detect_final(:,2),start_time,end_time,'VariableNames',var_names);
        
        [~,out_file,~]=fileparts(audiofile);
        writetable(out,strcat(out_dir,slash,out_file,'_detect_review.csv'));
        
        disp('Finish detection!');
    end
    
end
