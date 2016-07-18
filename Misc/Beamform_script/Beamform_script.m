%% Experiment Foldername %%
foldername ='/data/cfudata7/cavh/phantoms/microbubble/initial_try/SARUS/exp_2015.07.09_12.31';
my_sarus=sarus_type(foldername); % Load experiments parameters

log_compression = @(x) 20*log10(x/max(x(:)));
my_sarus = my_field;
%% Beamforming parameters  %%
% Axial parameters %
start_depth= 0.01; % Starting depth to beamform [meters]
end_depth= 0.03;   % Ending depth to beamform [meters]
lambda=my_sarus.phantom.c/my_sarus.xdc.fc;  % Wavelength in meters
dr = 50*10^(-6); % lambda/8;   % Resolution on the depth axis
% Lateral parameters %
start_lateral = -0.01; % my_sarus.xdc.element_positions(1,1);  % Start with leftmost element
end_lateral = 0.01; % my_sarus.xdc.element_positions(end,1);  % End with rightmost element
dr_lateral =  50*10^(-6);%(end_lateral-start_lateral)/1024;        % Resolution on the lateral dimension
% Time parameters %
seq_idx  = 1 ;  % Seq 1 is B-mode emission, seq 2 is a 5 emission faster flow sequence
frame_idx = 0;  % 1 fram have 128(1 b-mode) + 5*128(128 flow) emissions. Frame to beamform (zero = All)
em_idx  = 0;    % Emission to beamform (zero = All). Each picture from flow contains 5 emissions.


%% Cluster parameters %%
enable_cluster=true;
cluster_address= 'fcfu7';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Beamform      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% lines for BFT3 %
n_lines=length(start_lateral:dr_lateral:end_lateral); % Creates grid
lines.origin=zeros(n_lines,3);
lines.origin(:,1)=start_lateral:dr_lateral:end_lateral;
lines.origin(:,3)=lines.origin(:,3)+start_depth;
lines.origin=lines.origin';
lines.direction=repmat([0 0 1],[n_lines 1]); % Line direction
lines.direction=lines.direction';
lines.length=end_depth-start_depth; % Total span in depth
lines.dr=dr; % resolution in lateral direction

% Apodization for BFT3 %
xmt_apod.window='Tukey'; %xmt -> Transmit
xmt_apod.window_parameter=0.85;
xmt_apod.f=2;

rcv_apod.window='Hann'; % rcv -> Receive
rcv_apod.f=1;

xmt_apod.dynamic=true;
rcv_apod.dynamic=true;

if (enable_cluster)
    beamform_em(my_sarus,'seq_idx',seq_idx,...
        'frame_idx',frame_idx,...
        'em_idx', em_idx,...
        'name','BF_res1010_full_new',...
        'lines',lines,...
        'xmt_apod',xmt_apod,'rcv_apod',rcv_apod,...
        'hp_filter',300e3,...
        'do_save',true,...
        'cluster_node',cluster_address);        
else
    % Only one image is generated
    bft_image.nthreads=16;
    
    % Beamforming lines %
    rf_img=[];
    
    rf_line=beamform_em(my_sarus,'seq_idx',seq_idx,...
        'frame_idx',frame_idx,...
        'em_idx', em_idx,...
        'name','BF',...
        'lines',lines,...
        'xmt_apod',xmt_apod,'rcv_apod',rcv_apod,...
        'bft_image',bft_image,...
        'hp_filter',300e3,...
        'do_save',false);
    rf_img=rf_line{1}.bf_data;
    
    
    %Envelope detect and log compress
    compressed_data=log_compression(abs(rf_img));
    figure();
    %imagesc(compressed_data);
    %colormap(gray)
   
    
    
    imagesc([start_lateral end_lateral]*1000,[start_depth end_depth]*1000,compressed_data,[-dBmax 0]);
    colormap(gray)
    xlabel('Lateral distance x [mm]')
    ylabel('Depth z [mm]')
    axis('image')
    colorbar
end

%%

