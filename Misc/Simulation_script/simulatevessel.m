%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Simulation script for a vessel using flow sequence        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all;

fs  = 100e6;  % Field II Simulation fs


%% XDC & Sequence %%
% loads the 8670 and the default flow sequence
disp('Loading Sequence parameters.mat')

load('/home/cavh/Documents/MATLAB/data/sequences/bk8670_steered_defocused_chirp_wide/sequence');

%% Phantoms %%
% Vessel parameters
angle= 0;                 % Angle to simulate. 0 is axial 90 is transverse flow
no_frames= 10;           % Number of frames to simulate
center=[0.01 0 0.01];  % Center of the vessel 
fprf=1000;                % PRF to use
v_peak=-0.02;             % Peak velocity of the flow
c=1540;                   % Sound speed 
d_vessel=0.002;          % Diameter of the vessel 
length_multiplier=20;     % Length vessel multplier, length= multiplier* vessel diameter 
scat_density=1;        % Scatterer density [m^-3] (20e9 good for 7 MHz probe)
                          
                          
par.fname_base = '/data/cfudata6/s134082/Bachelorprojekt/simulation_data'; % Path to saving data 

my_sequence(1) = get_sequence_manual_delays(my_sequence(1),'xdc_pos',my_xdc.element_positions,'c',c);
my_sequence(2) = get_sequence_manual_delays(my_sequence(2),'xdc_pos',my_xdc.element_positions,'c',c);

my_sequence(1).prf=fprf;
my_sequence(2).prf=fprf;

  par.name = ['SA_flow_' num2str(angle) 'deg'];
    flow_angle = angle;
    % Create Phantom Object
    my_phantom = phantom_type(par);
    my_phantom.c=c;
    my_phantom=vessel_phantom(my_phantom,my_shuffle,...
        'center',center, ...
        'flow_angle',flow_angle,...
        'fprf',fprf,...
        'v_peak',v_peak,...
        'flow_type','parabolic',...
        'nframes',no_frames,...
        'd_vessel',d_vessel,...
        'sequence',2,...   
        'length_multiplier',length_multiplier,...
        'scat_density', scat_density);
    
    
    
    
    %%%%%   Create  FIELD II Object   %%%%%
    my_field = field_type('sample_freq',fs,...       
        'xdc', my_xdc, ...              % Transducer object 
        'phantom', my_phantom, ...      % Phantom object 
        'sequence',my_sequence,...      % Sequence object array
        'Frame_type','SA',... % Only stored
        'seq_shuffle', my_shuffle);     % The order the emissions in the sequence should be fired
    my_field.no_frames=no_frames;
    save_param(my_field);
    calc_scat_multi(my_field,'seq_idx',2);
    calc_scat_multi(my_field,'seq_idx',2,'merge','static');
    
    
    
    


