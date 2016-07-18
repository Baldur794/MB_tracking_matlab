%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Simulation script for a spinning disk                     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all;

fs  = 100e6;  % Field II Simulation fs
Seq_idx =  1;    % Sequence to simulate
no_frames = 5; % Number of frames to simulate
fprf=50;                % PRF to use

%% XDC & Sequence %%
% loads the 8670 and the default flow sequence
disp('Loading Sequence parameters.mat')

load('/home/cavh/Documents/MATLAB/data/sequences/bk8670_steered_defocused_chirp_wide/sequence');

%% Phantoms %%

%% Spinning disk phantom %%
center=[0 0 20]/1000;   % center [m]
v_peak=[0.0025];           % peak velocity [m/s]
radius=0.01;          % radius [m]
c= 1540;                % Sound speed [m/s]
                         
                          
par.fname_base = '/data/cfudata6/s134082/Bachelorprojekt/simulation_data/spinning_disk'; % Path to saving data 
par.name = ['Spinningdisk_' num2str(v_peak) 'm_s'];


% dumb_idx=find(my_shuffle(1,:)==Seq_idx);
% dumb_temp=my_shuffle(2,dumb_idx);
% dumb_temp= unique(dumb_temp); % for this shuffle
no_emissions=1*no_frames;

% Create Phantom Object
my_phantom = phantom_type(par);
my_phantom.c=c;

my_phantom=spinningdisk_phantom(my_phantom,...
    'center',center, ...
    'radius',radius,...
    'fprf',fprf,...
    'v_peak',v_peak,...
    'n_emissions',no_emissions...
    );


%my_sequence(1) = get_sequence_manual_delays(my_sequence(1),'xdc_pos',my_xdc.element_positions,'c',c);
my_sequence(1) = get_sequence_manual_delays(my_sequence(2),'xdc_pos',my_xdc.element_positions,'c',c);
my_sequence(2)=[];

my_sequence.v_sources(5)=[];
my_sequence.v_sources(4)=[];
my_sequence.v_sources(2)=[];
my_sequence.v_sources(1)=[];

my_shuffle = 3*ones(1,32);

my_sequence(1).prf=fprf;
%my_sequence(2).prf=fprf;

    %%%%%   Create  FIELD II Object   %%%%%
    my_field = field_type('sample_freq',fs,...       
        'xdc', my_xdc, ...              % Transducer object 
        'phantom', my_phantom, ...      % Phantom object 
        'sequence',my_sequence,...      % Sequence object array
        'Frame_type','SA');%,... % Only stored
        %'seq_shuffle', my_shuffle);     % The order the emissions in the sequence should be fired
my_field.no_frames=no_frames;
save_param(my_field);

%% Simulate %%
calc_scat_multi(my_field,'seq_idx',Seq_idx,'cluster_node','fcfu7');  

