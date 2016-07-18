function out = movement_tracking_b_mode_func( img_wind_cord, index )
%% Displacement RF
max_mov_y = 16;
max_mov_x = 4;

fps = 50;
frames = 250;
n_frames_rep = 50;
%start_frame = 103;
start_frame_new = 100;
H = block_matching;
H.max_mov_y = max_mov_y;
H.max_mov_x = max_mov_x;
H.cost_function = 'Norm Corr';
H.img_wind_cord = img_wind_cord;
idx_frame = 1;
mode = 'b_mode';

%% Displacement             
mov_y = [];
mov_x = [];
mov_y_list = zeros(n_frames_rep,frames);
mov_x_list = zeros(n_frames_rep,frames);
idx_mov = 1;
tempy = zeros(2,1);
tempx = zeros(2,1);

% Run through all templates
for idx_wind = 1:size(img_wind_cord,1)
    for start_frame = start_frame_new:start_frame_new+n_frames_rep+1
        % Load ref img
        idx_frame = start_frame;
        img_ref_wind = H.img_reference_window(idx_frame, idx_wind, mode);
        % Repeat for all frames
        for idx_frame = start_frame:start_frame+frames-1;
            %         % Load ref img
%             img_ref_wind = H.img_reference_window(idx_frame, idx_wind, mode);
            
            % Load new img
            img_new_temp = H.img_new_template(idx_frame, idx_wind, mode);
            
            [motion_y motion_x] = H.motion_displacement(img_ref_wind,img_new_temp);
            mov_y(idx_wind,idx_frame-start_frame+1) = motion_y;
            mov_x(idx_wind,idx_frame-start_frame+1) = motion_x;
        end
        
        if any(mov_y_list(1,:))
            [max_corr, lag] = crosscorr(mov_y,mov_y_list(1,:),60);
            lag_index = find(max_corr == max(max_corr(ceil(size(max_corr,2)/2):end)));
            tempy = [tempy [max(max_corr(ceil(size(max_corr,2)/2):end)); lag_index]];
            max_displacement = lag(lag_index)
            mov_y_list(idx_mov,max_displacement+1:end) = mov_y(1:end-max_displacement);
            
        else
            mov_y_list(1,:) = mov_y;
        end
        
        if any(mov_x_list(1,:))
            [max_corr, lag] = crosscorr(mov_x,mov_x_list(1,:),60);
            lag_index = find(max_corr == max(max_corr(ceil(size(max_corr,2)/2):end)));
            tempx = [tempx [max(max_corr(ceil(size(max_corr,2)/2):end)); lag_index]];
            max_displacement = lag(lag_index)
            mov_x_list(idx_mov,max_displacement+1:end) = mov_x(1:end-max_displacement);
            
        else
            mov_x_list(1,:) = mov_x;
        end
        
        idx_mov = idx_mov + 1;
    end
    mov_y = mean(mov_y_list);
    mov_y = mov_y(:,50:end);
    mov_x = mean(mov_x_list);
    mov_x = mov_x(:,50:end);
end
%% Calculate Veloctity (for fixed ref img)
% vel_x = [];
% vel_y = [];
% vel_abs = [];
% for idx_wind = 1:size(img_wind_cord,1)
%     vel_x(idx_wind,:) = mov_x(idx_wind,2:end)-mov_x(idx_wind,1:end-1);
%     vel_y(idx_wind,:) = mov_y(idx_wind,2:end)-mov_y(idx_wind,1:end-1);
%     vel_abs(idx_wind,:) = sqrt(vel_x(idx_wind,:).^2+vel_y(idx_wind,:).^2);
% end 

%% Calculate fft
% mov_y_frq = [];
% for idx_wind = 1:size(img_wind_cord,1)
%     mov_y_frq(idx_wind,:) = fft(mov_y(idx_wind,:)-mean(mov_y(idx_wind,:)));
% end
% 
% mov_x_frq = [];
% for idx_wind = 1:size(img_wind_cord,1)
%     mov_x_frq(idx_wind,:) = fft(mov_x(idx_wind,:)-mean(mov_x(idx_wind,:)));
% end

% %% Interpolate
% interpolate_factor = 11.51;
% mov_yy = spline(1:size(mov_y,2),mov_y,1:1/interpolate_factor:size(mov_y,2));
% mov_xx = spline(1:size(mov_x,2),mov_x,1:1/interpolate_factor:size(mov_x,2));
% 
% 
% % Axial displacement list
% temp = [];
% for h = 450:510
%     f_rep = 493;
%     axial_disp_list = zeros(idx_wind,f_rep,floor(size(mov_yy,2)/f_rep));
%     
%     for idx_wind = 1:size(img_wind_cord,1)
%         for i = 1:size(axial_disp_list,2)
%             for k = 1:size(axial_disp_list,3)
%                 index = size(axial_disp_list,2) * (k-1) + i;
%                 axial_disp_list(idx_wind,i,k) = mov_yy(idx_wind,index);
%             end
%         end
%     end
%     
%     
%     % Axial mean
%     axial_disp_mean = zeros(size(img_wind_cord,1),size(axial_disp_list,2));
%     axial_disp_mean = sum(axial_disp_list,3)/size(axial_disp_list,3);
%     
%     % Axial variance
%     axial_disp_var = zeros(size(img_wind_cord,1),1);
%     for idx_temp = 1:size(img_wind_cord,1)
%         axial_disp_var(idx_temp) = sum(var(axial_disp_list(idx_temp,:,:),1,3))/size(axial_disp_list,2);
%     end
%     %
%     % Axial std
%     axial_disp_std = zeros(size(img_wind_cord,1),1);
%     axial_disp_std = sqrt(axial_disp_var);
%     temp(h-449,:) = axial_disp_std ;
% end
% [t1 t2] = find(temp == min(temp(:)))
% 
% mov_y_mean = spline(linspace(1,f_rep/interpolate_factor,size(axial_disp_mean,2)),axial_disp_mean,1:f_rep/interpolate_factor);
% 
% 
% % Lateral displacement list
% temp = [];
% for h = 450:510
%     f_rep = 493;
%     lateral_disp_list = zeros(idx_wind,f_rep,floor(size(mov_xx,2)/f_rep));
%     
%     for idx_wind = 1:size(img_wind_cord,1)
%         for i = 1:size(lateral_disp_list,2)
%             for k = 1:size(lateral_disp_list,3)
%                 index = size(lateral_disp_list,2) * (k-1) + i;
%                 lateral_disp_list(idx_wind,i,k) = mov_xx(idx_wind,index);
%             end
%         end
%     end
%     
%     % lateral mean
%     lateral_disp_mean = zeros(size(img_wind_cord,1),size(lateral_disp_list,2));
%     lateral_disp_mean = sum(lateral_disp_list,3)/size(lateral_disp_list,3);
%     
%     % lateral variance
%     lateral_disp_var = zeros(size(img_wind_cord,1),1);
%     for idx_temp = 1:size(img_wind_cord,1)
%         lateral_disp_var(idx_temp) = sum(var(lateral_disp_list(idx_temp,:,:),1,3))/size(lateral_disp_list,2);
%     end
%     %
%     % lateral std
%     lateral_disp_std = zeros(size(img_wind_cord,1),1);
%     lateral_disp_std = sqrt(lateral_disp_var);
%     temp(h-449,:) = lateral_disp_std ;
% end
% [t3 t4] = find(temp == min(temp(:)))
% 
% mov_x_mean = spline(linspace(1,f_rep/interpolate_factor,size(lateral_disp_mean,2)),lateral_disp_mean,1:f_rep/interpolate_factor);
%% 
% locs = [];
% pks = [];
% for idx_wind = 1:size(img_wind_cord,1)
%     corr_match = xcorr(mov_y(idx_wind,:),mov_y_mean(idx_wind,:))
%     figure(); plot(corr_match);
%     findpeaks(corr_match,'NPeaks',floor(size(mov_y,2)/size(mov_y_mean,2)),'SortStr','descend','MinPeakDistance',10);
%     % Find peaks
%     [pks(idx_wind,:),locs(idx_wind,:)] = findpeaks(corr_match,'NPeaks',floor(size(mov_y,2)/size(mov_y_mean,2)),'SortStr','descend','MinPeakDistance',10);
% end
% locs = fliplr(locs)

% %%
% % Make compensation array axial
% mov_y_comp = zeros(size(img_wind_cord,1),size(mov_y,2));
% for idx_wind = 1:size(img_wind_cord,1)
%     for i = 1:size(locs,2)
%         mov_y_comp(idx_wind,locs(idx_wind,i)-size(mov_y,2)+1:locs(idx_wind,i)-size(mov_y,2)+size(mov_y_mean,2)) = mov_y_mean(idx_wind,:);
%     end
%     if size(mov_y_comp,2) > frames
%         mov_y_comp(:,frames+1:end) = [];
%     end
%     figure();
%     plot(mov_y_comp(idx_wind,:));
%     hold on
%     plot(mov_y(idx_wind,:));
%     hold off
%     xlabel('Frames'); ylabel('Axial displacement (mm)');
% end
% 
% % figure();
% % plot(abs(mov_y_comp(1,:)-mov_y(1,:)));
% 
% % Make compensation array lateral
% mov_x_comp = zeros(size(img_wind_cord,1),size(mov_x,2));
% for idx_wind = 1:size(img_wind_cord,1)
%     for i = 1:size(locs,2)
%         mov_x_comp(idx_wind,locs(idx_wind,i)-size(mov_x,2)+1:locs(idx_wind,i)-size(mov_x,2)+size(mov_x_mean,2)) = mov_x_mean(idx_wind,:);
%     end
%     if size(mov_x_comp,2) > frames
%         mov_x_comp(:,frames+1:end) = [];
%     end
%     figure();
%     plot(mov_x_comp(idx_wind,:));
%     hold on
%     plot(mov_x(idx_wind,:));
%     hold off
% end
out = {mov_y mov_x index};
end

