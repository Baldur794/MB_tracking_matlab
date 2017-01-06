%%%% Script for reading buffer dump from bk3000 %%%%
% The scripts reads the data, show the envelope and plot the RF center line
% of each stream

%log_compression = @(x) 20*log10(x/max(x(:)));
log_compression = @(x) 20*log10(x/5000);
addpath('./library/BkmDataFile');

folder_number=2;

%% Data location %%
% Contrast phantom %
%folderName='../../data/New_sequence_phantoms/MemDumps/10V-fc4.2M-bw1M';
%folderName='../../data/New_sequence_phantoms/MemDumps/40V-fc4.2M-bw1M';
%folderName='../../data/New_sequence_phantoms/MemDumps/40V-nofilter';
%folderName='../../data/New_sequence_phantoms/input-2016-11-15-110634-10V_rfdata';
%folderName='sequence_testing_phantom/input 2016-11-21-143447';
load('../../data/Sequence_testing_different_settings/folderNameList.mat'); % List with the seperate folder names
folderName = ['../../data/Sequence_testing_different_settings/' folderNameList{folder_number}];
folderNameList{folder_number}  % To display the folder
%folderName = '../../data/streaming_test/static 2min/stream streaming 2016-11-25-142218';

folderName = '/data/cfudata3/ramosh/cavh/microbubble-experiments/Part1-10V-6_5MHz';
%folderName = '/home/cavh/Documents/papers/Journals/SuperResolution/data/streaming_test/streaming 2016-12-05-143937';
%folderName = '/data/cfudata3/ramosh/cavh/microbubble-experiments/Part4-10V-4_5MHz-1uLs-1to100';
% Wire phantom %
% The phantom was not completely steady, misalignments occured.
%folderName='../../data/New_sequence_phantoms/input-2016-11-11-095056';
folderName ='/data/cfudata3/ramosh/cavh/microbubble-experiments/shottest/mod+--';
folderName ='/data/cfudata3/ramosh/cavh/microbubble-experiments/Part5-10V_6_5MHz_2Beam1Line_Flow_Ramp_130uls_to_0uls_60s_1to50';
%folderName = '/data/cfudata3/ramosh/cavh/microbubble-experiments/Part5-10V_6_5MHz_2Beam2Line_Flow_Ramp_130uls_to_0uls_60s_1to50';
%% Streams %%
streamName{1}='export-0.bkmdf';
streamName{2}='export-1.bkmdf';
streamName{3}='export-2.bkmdf';
% streamName{1}='stream-0.bkmdf';
% streamName{2}='stream-1.bkmdf';
% streamName{3}='stream-2.bkmdf';

sample_init= 100; % For removing the initial part affected by the excitation
sample_end= 1600; % For removing the bottom of the phantom
line_offset=0; % Line offset from the center line for SNR

no_frames=10;
no_skip_frames=120;
c = 1540;


%% Bandpass Filter parameters %%
bpf.filt_order=256;
bpf.cuttoff1= 2e6;            % Lower cutoff MHz
bpf.cuttoff2= 18e6;            % Higher cutoff MHz

%% Match Filter parameters
impulse_response_loc='/data/cfu/transducers/BK_9009/BK9009_impulse_response.mat'; % For match filtering

%% Read usecase %% (TODO)
[~, usecase, ext] = fileparts(folderName);
usecase = [usecase ext];
% sampling_rate= 60e6; % should be obtained from the usecase
% no_frames=15;
% c = 1540;
% start_line = 0;
% end_line = 149;
% element_pitch = 160*10^(-6);
% line_density = 1;
usecase_xml = xmlread([folderName '/usecase.xml']);
% read sampling rate
temp = usecase_xml.getElementsByTagName('receiveSampleFrequency');
sampling_rate = str2double(temp.item(0).getFirstChild.getData);
% read start line
temp = usecase_xml.getElementsByTagName('actualStartLine');
start_line = str2double(temp.item(0).getFirstChild.getData);
% read end line
temp = usecase_xml.getElementsByTagName('actualStopLine');
end_line = str2double(temp.item(0).getFirstChild.getData);
% read element pitch
temp = usecase_xml.getElementsByTagName('pitch');
element_pitch = str2double(temp.item(0).getFirstChild.getData);
% read line density
temp = usecase_xml.getElementsByTagName('lineDensity');
line_density = str2double(temp.item(0).getFirstChild.getData);
% read wave sample to calculate f0
temp = usecase_xml.getElementsByTagName('sampleCount');
f0_theory = (str2double(temp.item(0).getFirstChild.getData)/(120*10^6*1.5))^(-1);
temp = usecase_xml.getElementsByTagName('curveData');
waveform=getwaveform(temp.item(0).getFirstChild.getData);





% Load first frame from each stream
temp=[];
temp_data=[];
for it_stream=1:length(streamName)
    for it_frame=1:no_frames
        [header temp]=bkmload([folderName filesep streamName{it_stream}],270+it_frame);
        temp_data(:,:,it_frame) = hilbert(double(temp(sample_init:sample_end,:)));
        timestamp_frames(it_frame,it_stream)=header.timeStamp;
        lineCount_frames(it_frame,it_stream)=header.lineCount;
    end
    all_data{it_stream} = temp_data;
end

% Filter data, average and hilbert transform
%hpf.filter = designfilt('highpassfir', 'FilterOrder', 256, 'CutoffFrequency', 100000, 'SampleRate', 60000000);
%bpf.filter = designfilt('bandpassfir', 'StopbandFrequency1', 500000, 'PassbandFrequency1', 1000000, 'PassbandFrequency2', 12000000, 'StopbandFrequency2', 14000000, 'StopbandAttenuation1', 60, 'PassbandRipple', 1, 'StopbandAttenuation2', 60, 'SampleRate', 60000000);
bpf.filter = designfilt('bandpassfir', 'FilterOrder', bpf.filt_order, 'CutoffFrequency1', bpf.cuttoff1, 'CutoffFrequency2', bpf.cuttoff2, 'SampleRate', sampling_rate);

for it_stream=1:length(streamName)
    all_data{it_stream} = filtfilt(bpf.filter,all_data{it_stream});
    %     for idx_line = 1:size(all_data{it_stream},2)
    %         all_data{it_stream}(:,idx_line) = conv(all_data{it_stream}(:,idx_line),bpf.coef,'same');
    %     end
    
    data{it_stream}= mean(all_data{it_stream},3);
    data_SD{it_stream}= std(all_data{it_stream},0,3);
end

% Axis of data
t_axis=linspace(sample_init*(1/sampling_rate),sample_init*(1/sampling_rate)+(1/sampling_rate)*(size(data{1},1)-1),size(data{1},1)); % in seconds
fourier_axis=linspace(-sampling_rate/2,sampling_rate/2,size(data{1},1)); %in Hz
% start_axial = sample_init*(1/sampling_rate)*c*10^3; %in mm
%end_axial = sample_init*(1/sampling_rate)+(1/sampling_rate)*(size(data{1},1)-1)*c*10^3; %in mm
depth_axis = (t_axis*c)/2;   %  Values always in MKS
lateral_axis = ((start_line:end_line)-floor((end_line-start_line)/2))*element_pitch/line_density;
%start_lateral = 0; %in m
%end_lateral = (end_line-start_line)*element_pitch/line_density; % in [m]
% Number of ticks (spatial axis)
x_ticks = 4;
y_ticks = 5;

%%% Get the match filter  %%%
load(impulse_response_loc);
waveform(isnan(waveform))=[];
impulse_response=resample(imp.impulse_response,120e6,imp.fs); % TODO: load 120 from usecase
matched_filter= conv(conv(waveform,impulse_response,'same'),impulse_response); % Convolved twice
matched_filter=matched_filter(1:length(waveform));
matched_filter=matched_filter(end:-1:1);
matched_filter=resample(matched_filter,sampling_rate,120e6);


%Show the enveloped detected data
figure();
for it_streams=1:length(streamName)
    % Match filter 
    %rf_data=conv2(matched_filter,[1],data{it_streams},'same');
    rf_data=data{it_streams};
    subplot(1,length(streamName),it_streams)
    imagesc(lateral_axis*1000,depth_axis*1000,log_compression(abs(rf_data)),[-60 0])
    title([usecase ' ' streamName{it_streams}])
    xlabel('lateral [mm]')
    axis image
    % Is better to spply the axis directly
    %     set(gca,'XTick',linspace(1,size(data{1},2),x_ticks))
    %     set(gca,'XTickLabel',round(linspace(start_lateral,end_lateral,x_ticks)))
    ylabel('depth [mm]')
    %     set(gca,'YTick',linspace(1,size(data{1},1),y_ticks))
    %     set(gca,'YTickLabel',round(linspace(start_axial,end_axial,y_ticks)))
end
colormap(gray(255))

% Plot the RFdata
 figure();
rf_data=[];
for it_streams=1:length(streamName)
    rf_data(:,it_streams)=data{it_streams}(:,floor(size(data{it_streams},2)/2)+line_offset);
end

subplot(2,3,[1 2])
plot(t_axis*1e6,real(rf_data))
xlim([t_axis(1) t_axis(end)]*1e6)
xlabel('time [us]')
title([usecase ' RF single line all streams'])
subplot(2,3,3)
plot(fourier_axis/1e6,abs(fftshift(fft(real(rf_data),[],1))))
xlim([0 20])
title(['Fourier RF single line all streams'])


subplot(2,3,[4 5])
plot(t_axis*1e6,real(rf_data(:,1))+real(rf_data(:,2)),t_axis*1e6,real(rf_data(:,2))-real(rf_data(:,3)))
title('PI sequence single line')
xlabel('time [us]')
xlim([t_axis(1) t_axis(end)]*1e6)
subplot(2,3,[6])
plot(fourier_axis/1e6,abs(fftshift(fft(real(rf_data(:,1))+real(rf_data(:,2)),[],1))), ...
    fourier_axis/1e6,abs(fftshift(fft(real(rf_data(:,3))+real(rf_data(:,2)),[],1))))
xlim([0 20])
title(['Fourier PI single line all streams'])


% Show envelope of added streams
for it_stream=1:length(streamName)-1
    PI_data{it_stream}=mean(all_data{it_stream}+all_data{it_stream+1},3);
    PI_data_SD{it_stream}=std(all_data{it_stream}+all_data{it_stream+1},[],3);
end

figure();
load(impulse_response_loc);
waveform(isnan(waveform))=[];
impulse_response=resample(imp.impulse_response,120e6,imp.fs); % TODO: load 120 from usecase
matched_filter= conv(conv(waveform,impulse_response,'same'),impulse_response); % Convolved twice
matched_filter=matched_filter(1:length(waveform));
matched_filter=matched_filter(end:-1:1);
matched_filter=resample(matched_filter,sampling_rate,120e6);

for it_stream=1:length(streamName)-1
    subplot(1,length(streamName)-1,it_stream)
    %rf_data=conv2(matched_filter,[1],PI_data{it_stream},'same');
    rf_data=PI_data{it_stream};
    compressed_PI=log_compression(abs(rf_data));
    imagesc(lateral_axis*1000,depth_axis*1000,compressed_PI,[-30 0])
    title([usecase ' PI sequence image'])
    xlabel('lateral [mm]')
    axis image
    %     set(gca,'XTick',linspace(1,size(data{1},2),x_ticks))
    %     set(gca,'XTickLabel',round(linspace(start_lateral,end_lateral,x_ticks)))
    ylabel('depth [mm]')
    %     set(gca,'YTick',linspace(1,size(data{1},1),y_ticks))
    %     set(gca,'YTickLabel',round(linspace(start_axial,end_axial,y_ticks)))
end
colormap(gray(255))


% Plot SNR
% Average on a selected region(middel)
win_len  = 100;
lpf= fir1(win_len,0.01)';
lines_to_avg=line_offset+(floor(size(data{it_streams},2)/2)-15:1:floor(size(data{it_streams},2)/2)+15);
for it_stream=1:length(streamName)
    SNR(:,it_stream) = 20*log10(mean(abs(data{it_stream}(:,lines_to_avg)./data_SD{it_stream}(:,lines_to_avg)),2));
    SNR(:,it_stream) = conv(SNR(:,it_stream), lpf, 'same');
end

for it_stream=1:length(streamName)-1
    SNR_PI(:,it_stream) = 20*log10(mean(abs(PI_data{it_stream}(:,lines_to_avg)./PI_data_SD{it_stream}(:,lines_to_avg)),2));
    SNR_PI(:,it_stream) = conv(SNR_PI(:,it_stream), lpf, 'same');
end

% 
% figure();
% subplot(2,1,1)
% plot((t_axis*c)/2*1000,SNR,'linewidth',2)
% xlabel('depth [mm]')
% title('SNR Bmode signal')
% subplot(2,1,2)
% plot((t_axis*c)/2*1000,SNR_PI,'linewidth',2)
% xlabel('depth [mm]')
% title('SNR PI signal')

% % Display SNR through time
% display_snr=input('Display SNR through time?  [y] ' ,'s');
% if isempty(display_snr) %|| strcomp(display_snr)
%     
%     for it_stream_window=1:20
%         
%         % Load the stream
%         frame_ini=(no_skip_frames*(it_stream_window-1))+1
%         temp=[];
%         temp_data=[];
%         for it_stream=1:length(streamName)
%             for it_frame=1:no_frames
%                 try
%                     [header temp]=bkmload([folderName filesep streamName{it_stream}],frame_ini+it_frame);
%                 catch
%                     break;
%                 end
%                 temp_data(:,:,it_frame) = hilbert(double(temp(sample_init:sample_end,:)));
%                 timestamp_frames(it_frame,it_stream)=header.timeStamp;
%             end
%             all_data{it_stream} = temp_data;
%         end
%         
%         % Filter data, average and hilbert transform
%         %hpf.filter = designfilt('highpassfir', 'FilterOrder', 256, 'CutoffFrequency', 100000, 'SampleRate', 60000000);
%         %bpf.filter = designfilt('bandpassfir', 'StopbandFrequency1', 500000, 'PassbandFrequency1', 1000000, 'PassbandFrequency2', 12000000, 'StopbandFrequency2', 14000000, 'StopbandAttenuation1', 60, 'PassbandRipple', 1, 'StopbandAttenuation2', 60, 'SampleRate', 60000000);
%         bpf.filter = designfilt('bandpassfir', 'FilterOrder', bpf.filt_order, 'CutoffFrequency1', bpf.cuttoff1, 'CutoffFrequency2', bpf.cuttoff2, 'SampleRate', sampling_rate);
%         for it_stream=1:length(streamName)
%             all_data{it_stream} = filtfilt(bpf.filter,all_data{it_stream});
%             %     for idx_line = 1:size(all_data{it_stream},2)
%             %         all_data{it_stream}(:,idx_line) = conv(all_data{it_stream}(:,idx_line),bpf.coef,'same');
%             %     end
%             data{it_stream}= mean(all_data{it_stream},3);
%             data_SD{it_stream}= std(all_data{it_stream},0,3);
%         end
%         % Show envelope of added streams
%         for it_stream=1:length(streamName)-1
%             PI_data{it_stream}=mean(all_data{it_stream}+all_data{it_stream+1},3);
%             PI_data_SD{it_stream}=std(all_data{it_stream}+all_data{it_stream+1},[],3);
%         end
%         
%         
%         %Estimate SNR curve
%         % Average on a selected region(middel)
%         win_len  = 100;
%         lpf= fir1(win_len,0.01)';
%         lines_to_avg=line_offset+(floor(size(data{it_streams},2)/2)-15:1:floor(size(data{it_streams},2)/2)+15);
%         it_stream=1;
%         SNR_temp_b=[];
%         SNR_temp_b = 20*log10(mean(abs(data{it_stream}(:,lines_to_avg)./data_SD{it_stream}(:,lines_to_avg)),2));
%         SNR_temp_b  = conv(SNR_temp_b, lpf, 'same');
%         
%         
%         it_stream=2;
%         SNR_temp_pi=[];
%         SNR_temp_pi = 20*log10(mean(abs(PI_data{it_stream}(:,lines_to_avg)./PI_data_SD{it_stream}(:,lines_to_avg)),2));
%         SNR_temp_pi = conv(SNR_temp_pi, lpf, 'same');
%         
%         
%         % Save SNRs
%         SNR_bmode_all(:,it_stream_window)= SNR_temp_b ;
%         SNR_PI_all(:,it_stream_window)= SNR_temp_pi ;
%         
%     end
%     
% end

%figure(); imagesc(1:20,(t_axis*c)/2*1000, SNR_PI_all,[-10 10])





