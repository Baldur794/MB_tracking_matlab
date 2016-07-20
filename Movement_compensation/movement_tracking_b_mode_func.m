function out = movement_tracking_b_mode_func( img_wind_cord, index )
%% Displacement RF

max_mov_y = 10;
max_mov_x = 10;

fps = 50;
frames = 1000;
n_frames_rep = 50;
start_frame = 1;
H = block_matching;
H.max_mov_y = max_mov_y;
H.max_mov_x = max_mov_x;
H.cost_function = 'xCorr';
H.img_wind_cord = img_wind_cord;
idx_frame = 1;
mode = 'b_mode';
             
%% Displacement             
mov_y = [];
mov_x = [];

% Run through all templates
for idx_wind = 1:size(img_wind_cord,1)
        % Load ref img
        idx_frame = start_frame;
        img_ref_wind = H.img_reference_window(idx_frame, idx_wind, mode);
        % Repeat for all frames
        for idx_frame = start_frame:start_frame+frames-1;
            %         % Load ref img
            img_ref_wind = H.img_reference_window(idx_frame, idx_wind, mode);
            
            % Load new img
            img_new_temp = H.img_new_template(idx_frame+1, idx_wind, mode);
            
            [motion_y motion_x] = H.motion_displacement(img_ref_wind,img_new_temp);
            mov_y(idx_wind,idx_frame-start_frame+1) = motion_y;
            mov_x(idx_wind,idx_frame-start_frame+1) = motion_x;
        end        
end
out = {mov_y mov_x index};
end

