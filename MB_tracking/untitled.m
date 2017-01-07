%%%% Script for reading streams from bk3000 %%%%
% The scripts reads the bk stream and makes a  %
% backgrounf filtered video of the PI sequence.%
%                                              %
% The scripts is meant to be tested on a       %
% flow phantom.                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'defaultaxesfontsize',14);
set(0,'defaulttextfontsize',14);

log_compression = @(x) 20*log10(x/5000);
% addpath('./library/BkmDataFile');


%% Data location %%
% flow phantom %
folderName ='/data/cfudata3/ramosh/cavh/microbubble-experiments/twochannel_5ul_sec';%%_sec_freehand';
% video_frames_folder=[folderName '/video/imgs'];
% Lost synch flag 
fix_synch=false;     % It seems that streams lost synch when using the full aperture, to
                    % correct enable this flag. 
                    
                       
                    
%% Streams %%
streamName{1}='export-0.bkmdf';
streamName{2}='export-1.bkmdf';

sample_init= 700; % For removing the initial part affected by the excitation  700-1000 Freehand
sample_end= 1000; % For removing the bottom of the phantom

no_frames=10;           % Number of frames for moving background removal
no_skip_frames_ini=10000;  % Number of Frames to skip at the beginning of the video
skip_frames=1;           % Number of frames to skip 
no_total_frames=10000;     % Number of frames in the video
c = 1540;                % sound speed [m/s]
               
%% Bandpass Filter parameters %%
bpf.filt_order=256;
bpf.cuttoff1= 2e6;            % Lower cutoff MHz
bpf.cuttoff2= 18e6;            % Higher cutoff MHz


%% Read usecase %% 
[~, usecase, ext] = fileparts(folderName);
usecase = [usecase ext];
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



%% Load DATA 
for it_vid_frame=1:no_total_frames
% Load first frame from each stream
temp=[];
temp_data=[];
timestamp_frames=[];
lineCount_frames=[];
all_data=[];
% Hack for synching 
if fix_synch
    no_frames=no_frames+1;  % Load an extra frame
end
    
for it_stream=1:length(streamName)
    for it_frame=1:no_frames
        [header temp]=bkmload([folderName filesep streamName{it_stream}],no_skip_frames_ini+(skip_frames*(it_vid_frame-1))+it_frame);
        temp_data(:,:,it_frame) = hilbert(double(temp(sample_init:sample_end,:)));
        timestamp_frames(it_frame,it_stream)=header.timeStamp;
        lineCount_frames(it_frame,it_stream)=header.lineCount;
    end
    all_data{it_stream} = temp_data;
end

% HACK 
if fix_synch
    all_data{1}(:,:,end)=[];  % Remove last one
    all_data{2}(:,:,1)=[];    % Remove first one
    timestamp_frames(:,2)= circshift(timestamp_frames(:,2),-1,1);
    timestamp_frames(end,:)=[];
    lineCount_frames(:,2)= circshift(lineCount_frames(:,2),-1,1);
    lineCount_frames(end,:)=[];
    no_frames=no_frames-1;  
end

%% Check streams consistency
% Check for dropped lines (frames) between streams 
temp_line_count=abs(diff(lineCount_frames,[],2))-1;
err_frames=find(temp_line_count~=0);
dropped_line_count= length(err_frames);

% Drop the erroneous frames
for it =1:dropped_line_count
    for it_stream=1:length(streamName)
        all_data{it_stream}(:,:,err_frames(end-it+1))=[];
    end
end
    


temp= all_data{1}+all_data{2} -repmat(mean(all_data{1}+all_data{2},3),[1 1 size(all_data{2},3)]);

vid_frame(:,:, it_vid_frame)= temp(:,:,floor(size(temp,3)/2));
end

% % Axis of data
% t_axis=linspace(sample_init*(1/sampling_rate),sample_init*(1/sampling_rate)+(1/sampling_rate)*(size(all_data{1},1)-1),size(all_data{1},1)); % in seconds
% depth_axis = (t_axis*c)/2;   %  Values always in MKS
% lateral_axis = ((start_line:end_line)-floor((end_line-start_line)/2))*element_pitch/line_density;
% fighandle=figure();
% for it=1:9146
figure(2)
    imagesc(lateral_axis*1000,depth_axis*1000,abs(vid_frame(:,:,1)),[0 2000])
    colormap(gray(255))
%     axis image
%     xlabel('Lateral [mm]')
%     ylabel('depth [mm]')
%     name = sprintf([ video_frames_folder '/frame%04d'] , it);
%     print(fighandle, '-dpng', '-r300', name)
% end
