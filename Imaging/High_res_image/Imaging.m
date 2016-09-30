%% Age and Density filter
MB_index_filter_frame = [];
MB_index_filter_age = [];
MB_index_filter_count = [];
MB_index_filter_stuck_vel = []; % Filtered list of MB's satisfying conditions
MB_index_filter_stuck_int = [];
MB_index_filter_avg = []; % Filtered list of MB's satisfying conditions
MB_index_filter_dir = []; % Filtered list of MB's satisfying conditions

MB_index_list = []; % List of indexes with for found MB's
MB_vel_list = []; % List with all velocity vectors
temp = [];
MB_log_copy = MB_log;
MB_index_init = 1:size(MB_log_copy,2);
% 
% MB_log_copy = [];
% MB_log_copy = MB_plot(1);
% MB_log_copy(2) = MB_plot(3);
% MB_index_init = 1:size(MB_log_copy,2);

% check age condition
if MB_frame_condition_min ~= 0; 
    for i = 1:size(MB_index_init,2)
        MB_index = MB_index_init(i);
        if (MB_log_copy(MB_index).age(1) >= MB_frame_condition_min) && (MB_log_copy(MB_index).age(2) <= MB_frame_condition_max)
            MB_index_filter_frame = [MB_index_filter_frame, MB_index];
        end
    end
else
    MB_index_filter_frame = MB_index_init;
end

% check age condition
for i = 1:size(MB_index_filter_frame,2)
    MB_index = MB_index_filter_frame(i);
   if (MB_log_copy(MB_index).age(3) >= MB_age_condition_min) && (MB_log_copy(MB_index).age(3) <= MB_age_condition_max)
    MB_index_filter_age = [MB_index_filter_age, MB_index];
   end
end

% check count condition
for i = 1:size(MB_index_filter_age,2)
    MB_index = MB_index_filter_age(i);
   if (max(MB_log_copy(MB_index).count(:)) <= MB_count_condition)
    MB_index_filter_count = [MB_index_filter_count, MB_index];
   end
end

%------ 
%Remove stuck MB's from average speed
if MB_avg_vel_condition ~= 0
    for i = 1:size(MB_index_filter_count,2)
        MB_index = MB_index_filter_count(i);
        
        avg_vel = sqrt((MB_log_copy(MB_index).centroid(1,1)-MB_log_copy(MB_index).centroid(end,1))^2+(MB_log_copy(MB_index).centroid(1,2)-MB_log_copy(MB_index).centroid(end,2))^2)/MB_log_copy(MB_index).age(3);
        
        % Check density condition
        if (avg_vel >= MB_avg_vel_condition)%
            MB_index_filter_stuck_vel = [MB_index_filter_stuck_vel, MB_index];
        end
    end
else
    MB_index_filter_stuck_vel = MB_index_filter_count;
end
%------

MB_index_list = [];
% List all MB positions and their velocities
for i = 1:size(MB_index_filter_stuck_vel ,2)
   MB_index = MB_index_filter_stuck_vel(i);
   MB_index_list = [MB_index_list, sub2ind([img_size(1:2)],MB_log_copy(MB_index).centroid(1:end-1,1) ,MB_log_copy(MB_index).centroid(1:end-1,2))'];
end

% Scatter image from all MB locations
scatter_matrix = zeros([img_size(1:2)]); 
for i = 1:size(MB_index_list,2)
    scatter_matrix(MB_index_list(i)) = scatter_matrix(MB_index_list(i)) + 1;
end


% Find avg density of other MB for each MB location
for i = 1:size(MB_index_filter_stuck_vel ,2)
    MB_index = MB_index_filter_stuck_vel(i);
    MB_density_list = [];
    
    index_list_MB_index = [sub2ind([img_size(1:2)],MB_log_copy(MB_index).centroid(1:end-1,1) ,MB_log_copy(MB_index).centroid(1:end-1,2))']; % list of positions and vel for given MB
    
    % Scatter image without the current MB
    logical_matrix_MB_index = scatter_matrix;
    for j = 1:size(index_list_MB_index,2)
        logical_matrix_MB_index(index_list_MB_index(j)) = logical_matrix_MB_index(index_list_MB_index(j)) - 1; % Removes the MB in interest
    end
    
    % Calculated density of MB in window
    for j = 1:size(index_list_MB_index,2)
        % y,x subscripts from direct matrix indexes
        [subscript_y, subscript_x] = ind2sub([img_size(1:2)],index_list_MB_index(j));
        
        % Update window coordinates
        MB_window_coord(1,1) = subscript_y-MB_window_size_density_avg(1);
        MB_window_coord(2,1) = subscript_y+MB_window_size_density_avg(1);
        MB_window_coord(1,2) = subscript_x-MB_window_size_density_avg(2);
        MB_window_coord(2,2) = subscript_x+MB_window_size_density_avg(2);
        
        % Check for out of bounce
        % y1
        if MB_window_coord(1,1) <= 0
            MB_window_coord(1,1) = 1;
        end
        % y2
        if MB_window_coord(2,1) > img_size(1)
            MB_window_coord(2,1) = img_size(1);
        end
        % x1
        if MB_window_coord(1,2) <= 0
            MB_window_coord(1,2) = 1;
        end
        % x2
        if MB_window_coord(2,2) > img_size(2)
            MB_window_coord(2,2) = img_size(2);
        end
        
        % Sum number of MB in window
        logical_matrix_MB_index_window = logical_matrix_MB_index(MB_window_coord(1,1):MB_window_coord(2,1),MB_window_coord(1,2):MB_window_coord(2,2));
        MB_density_list(j) = sum(sum(logical_matrix_MB_index_window)); % number of MB's found around given MB with # MB_index
    end
    % Calculate average
    MB_density = sum(MB_density_list)/MB_log_copy(MB_index).age(3);
    
    % Check density condition
    if (MB_density >= MB_dens_condition_avg)% && (min(MB_density_list) >= MB_index_filter_dir)
        MB_index_filter_avg = [MB_index_filter_avg, MB_index];
    end
end

% Find dir
if kidney_center(1) ~= 0
    for i = 1:size(MB_index_filter_avg ,2)
        MB_index = MB_index_filter_avg(i);
        
        distance_to_center_start = sqrt((MB_log_copy(MB_index).centroid(1,1)-kidney_center(1))^2+(MB_log_copy(MB_index).centroid(1,2)-kidney_center(2))^2);
        distance_to_center_end = sqrt((MB_log_copy(MB_index).centroid(end,1)-kidney_center(1))^2+(MB_log_copy(MB_index).centroid(end,2)-kidney_center(2))^2);
        % Check density condition
        if (distance_to_center_start < distance_to_center_end)
            MB_index_filter_dir = [MB_index_filter_dir, MB_index];
        end
    end
else
    MB_index_filter_dir = MB_index_filter_avg;
end

 %----------------------------
MB_index_list = []; % List of indexes with for found MB's
% List all MB positions and their velocities
for i = 1:size(MB_index_filter_dir ,2)
    MB_index = MB_index_filter_dir(i);
    MB_index_list = [MB_index_list, sub2ind([img_size(1:2)],MB_log_copy(MB_index).centroid(1:end-1,1) ,MB_log_copy(MB_index).centroid(1:end-1,2))'];
end

% Scatter image from all MB locations
scatter_matrix = zeros([img_size(1:2)]); 
for i = 1:size(MB_index_list,2)
    scatter_matrix(MB_index_list(i)) = scatter_matrix(MB_index_list(i)) + 1;
end

 %-----------------------------
 n_vel_avg = 50;
 
 if n_vel_avg > MB_age_condition_min-2
     n_vel_avg = MB_age_condition_min-2;
 end
 % Smooth velocities
 for i = 1:size(MB_index_filter_single,2)
     MB_index = MB_index_filter_single(i);
    for j = 1:MB_log_copy(MB_index).age(3)-1
        if floor(j-(n_vel_avg-1)/2) < 1
            MB_log_copy(MB_index).vel(j,1) = sum(MB_log_copy(MB_index).vel(1:1+(n_vel_avg-1),1))/n_vel_avg;
            MB_log_copy(MB_index).vel(j,2) = sum(MB_log_copy(MB_index).vel(1:1+(n_vel_avg-1),2))/n_vel_avg;
        elseif floor(j+(n_vel_avg-1)/2) > MB_log_copy(MB_index).age(3)-1
            MB_log_copy(MB_index).vel(j,1) = sum(MB_log_copy(MB_index).vel(MB_log_copy(MB_index).age(3)-1-(n_vel_avg-1):MB_log_copy(MB_index).age(3)-1,1))/n_vel_avg;
            MB_log_copy(MB_index).vel(j,2) = sum(MB_log_copy(MB_index).vel(MB_log_copy(MB_index).age(3)-1-(n_vel_avg-1):MB_log_copy(MB_index).age(3)-1,2))/n_vel_avg;
        else
            MB_log_copy(MB_index).vel(j,1) = sum(MB_log_copy(MB_index).vel(floor(j-(n_vel_avg-1)/2):floor(j+(n_vel_avg-1)/2),1))/n_vel_avg;
            MB_log_copy(MB_index).vel(j,2) = sum(MB_log_copy(MB_index).vel(floor(j-(n_vel_avg-1)/2):floor(j+(n_vel_avg-1)/2),2))/n_vel_avg;
        end
    end
 end
 
% List of MB positions and their velocities satisfying conditions
MB_index_list = [];
MB_coord_list = [];
MB_vel_list = [];
for i = 1:size(MB_index_filter_dir,2)
  MB_index = MB_index_filter_dir(i);
  MB_index_list = [MB_index_list, sub2ind([img_size(1:2)],MB_log_copy(MB_index).centroid(2:end-1,1) ,MB_log_copy(MB_index).centroid(2:end-1,2))'];
  MB_coord_list = [MB_coord_list, MB_log_copy(MB_index).centroid(2:end-1,:)'];
  MB_vel_list = [MB_vel_list, [MB_log_copy(MB_index).vel(2:end,1),MB_log_copy(MB_index).vel(2:end,2)]'];
end
% MB_vel_list_abs = sqrt(MB_vel_list(1,:).^2+MB_vel_list(2,:).^2)* 10*10^(-6)*fps*10^3;
% MB_vel_list_dir = mod(atan2(MB_vel_list(1,:),MB_vel_list(2,:))+5/2*pi,2*pi);

% Scatter image from all MB locations
scatter_matrix = zeros([img_size(1:2)]); 
for i = 1:size(MB_index_list,2)
    scatter_matrix(MB_index_list(i)) = scatter_matrix(MB_index_list(i)) + 1;
end

% % Velocity img, Amplitude and direction
% % Final velocity and direction img
% vel_abs_img = zeros([img_size(1:2)]);
% vel_dir_img = zeros([img_size(1:2)]);
% 
% % Weighing filter for weighted velocities
% %weighing_filter_radius = 3; % Radius around centroid to be considered
% weighing_filter_size = 2*weighing_filter_radius+1; % Must be an odd number
% weighing_filter_halfsize = floor(weighing_filter_size/2); % Only to make indexing easier
% weighing_filter = zeros(weighing_filter_size); % Contains the final filter
% weighing_center = [ceil(weighing_filter_size/2),ceil(weighing_filter_size/2)]; % Center coordinates
% %weighing_factor = 1; % Distance weighing factor
% 
% % Makes linear weighing filter for full square
% for i = 1:size(weighing_filter,1)
%     for j = 1:size(weighing_filter,2)
%         weighing_filter(i,j) = sqrt((i-weighing_center(1))^2+(j-weighing_center(2))^2);
%     end
% end
% 
% % Discard indexes with larger radius than weighing_filter_radius
% if weighing_factor == 0
%     weighing_filter = zeros(size(weighing_filter));
%     weighing_filter(ceil(size(weighing_filter,1)/2),ceil(size(weighing_filter,2)/2)) = 1;
% else
%     temp_weighing_filter = weighing_filter;
%     temp_weighing_filter(weighing_filter > weighing_filter_radius) = 0;
%     temp_weighing_filter(weighing_filter <= weighing_filter_radius) = 1;
%     weighing_filter = weighing_filter * 1/weighing_factor;
%     weighing_filter = exp(-(weighing_filter/weighing_filter_radius).^2);
%     weighing_filter = weighing_filter.*temp_weighing_filter;
% end
% 
% % Velocity matrix
% vel_logical_matrix = zeros([img_size(1:2)]);
% vel_x = zeros([img_size(1:2)]); % x-part of velocity
% vel_y = zeros([img_size(1:2)]); % y-part of velocity
% 
% for i = 1:size(MB_index_list,2)
%     vel_x(MB_index_list(i)) = vel_x(MB_index_list(i)) + MB_vel_list(2,i);
%     vel_y(MB_index_list(i)) = vel_y(MB_index_list(i)) + MB_vel_list(1,i);
%     vel_logical_matrix(MB_index_list(i)) = vel_logical_matrix(MB_index_list(i)) + 1;
% end
% 
% % Weight multiple occuring positions
% vel_x(find(vel_x > 0)) = vel_x(find(vel_x > 0))./vel_logical_matrix(find(vel_x > 0));
% vel_y(find(vel_y > 0)) = vel_y(find(vel_y > 0))./vel_logical_matrix(find(vel_y > 0));
% vel_logical_matrix(vel_logical_matrix > 1) = 1;
% 
% % y,x subscripts from direct matrix indexes
% [subscript_y, subscript_x] = ind2sub([img_size(1:2)],MB_index_list);
% 
% % Temporary window and coordinates area in interst
% vel_window_coord = [];
% vel_window = []; 
% vel_window_y = [];
% vel_window_x = [];
% 
% for i = 1:size(MB_index_list,2)
%     
%     % Update window coordinates
%     vel_window_coord(1,1) = subscript_y(i)-weighing_filter_halfsize;
%     vel_window_coord(2,1) = subscript_y(i)+weighing_filter_halfsize;
%     vel_window_coord(1,2) = subscript_x(i)-weighing_filter_halfsize;
%     vel_window_coord(2,2) = subscript_x(i)+weighing_filter_halfsize;
%     
%     % Check for out of bounce
%     if (vel_window_coord(1,1) > 0) && (vel_window_coord(2,1) < img_size(1)) && (vel_window_coord(1,2) > 0) && (vel_window_coord(2,2) < img_size(2))
%         % Windows of interest around specific position
%         vel_window = vel_logical_matrix(vel_window_coord(1,1):vel_window_coord(2,1),vel_window_coord(1,2):vel_window_coord(2,2));
%         vel_window_y = vel_y(vel_window_coord(1,1):vel_window_coord(2,1),vel_window_coord(1,2):vel_window_coord(2,2));
%         vel_window_x = vel_x(vel_window_coord(1,1):vel_window_coord(2,1),vel_window_coord(1,2):vel_window_coord(2,2));
%         
%         % Distance weighting exp(-(d_i/r)^2)
%         vel_weight_matrix = weighing_filter.*vel_window;
%         
%         % Sum of all weightings
%         vel_weight_Z = sum(vel_weight_matrix(:)); 
%         
%         % Final weighted velocities
%         vel_weight_y = sum(sum(1/vel_weight_Z*vel_weight_matrix.*vel_window_y));
%         vel_weight_x = sum(sum(1/vel_weight_Z*vel_weight_matrix.*vel_window_x));
%         
%         % insert weighted velocities to right positions
%         vel_y(MB_index_list(i)) = vel_weight_y;
%         vel_x(MB_index_list(i)) = vel_weight_x;
%     end
% end
% 
% % Calculate absolute velocity
% vel_abs_img = sqrt(vel_y.^2+vel_x.^2);
% % Correct for pixel resolution 
% vel_abs_img = vel_abs_img * 10*10^(-6)*fps*10^3;
% 
% % Color map for direction img
% colormap_vel_dir = hsv;
% colormap_vel_dir(1,:) = [0 0 0];
% 
% % Direction img. Outputs in radians
% vel_dir_img = atan2(vel_y,vel_x);
% % From [-pi,pi] -> [0,2*pi]
% vel_dir_img(vel_dir_img ~= 0) = vel_dir_img(vel_dir_img ~= 0) + pi;
% % Turn colordisk
% color_rotate = 3/2*pi;
% vel_dir_img(vel_dir_img ~= 0) = mod(vel_dir_img(vel_dir_img ~= 0) + color_rotate,2*pi);
% 
% % Subscripts and index of non zero data pixels (both for vel and dir)
% [vel_abs_img_list_y vel_abs_img_list_x] = ind2sub(size(vel_abs_img),find(vel_abs_img(:) > 0)); 
% vel_abs_img_list = find(vel_abs_img(:));
% 
% 
% % List dir values at non-zero indexes
% vel_dir_img_list_val = vel_dir_img(vel_abs_img_list);
% % Adjust dir values to lay between [1:64] for colormap
% vel_dir_img_list_val = floor(vel_dir_img_list_val/(2*pi)*63+1);
% %----


%% Make direction image
% Make figure handle
fh = figure; clf;
% set position on screen
set(fh,'position',[-1850 570 560 420]);
subplot(1,2,1);

% Dir scatter plot
% colormapping for direction img
dir_color = hsv; 

% Show scatterplot with non-zero positions
fig = scatter(0, 0);
% fig = scatter(vel_abs_img_list_x, vel_abs_img_list_y);
% fig.CData = dir_color(vel_dir_img_list_val,:); % setting colors of individual dots depending on direction
MB_XData = [];
MB_YData = [];
MB_CData = [];

i_end = 1;
for i = 1:size(MB_index_filter_dir,2)
    MB_index = MB_index_filter_dir(i);
    
    i_start = i_end;
    i_end = i_end+MB_log_copy(MB_index).age(3)-3;
    
    MB_XData = [MB_XData MB_log_copy(MB_index).centroid(2:end-1,2)'];%scatter(vel_abs_img_list_x(i), vel_abs_img_list_y(i));
    MB_YData = [MB_YData MB_log_copy(MB_index).centroid(2:end-1,1)'];
    MB_CData = [MB_CData; dir_color(floor(mod(atan2(MB_log_copy(MB_index).vel(2:end,1),MB_log_copy(MB_index).vel(2:end,2))+5/2*pi,2*pi)/(2*pi)*63+1),:)];% setting colors of individual dots depending on direction
end
MB_index_list = [];
% MB_coord_list = [];
% MB_vel_list = [];
fig.XData = MB_XData;
fig.YData = MB_YData;
fig.CData = MB_CData;

% Scatterplot settings
set(fig,'SizeData', 0.5); % size of dots
set(fig,'MarkerFacecolor','flat'); % appearance of dots
xlabel('Lateral [mm]'); ylabel('Axial [mm]'); % title('Micro-Bubble image');
set(gca,'Xtick',linspace(50,1000,6)); set(gca, 'XTickLabel',linspace(0,10,6));
set(gca,'Ytick',linspace(60,1950,11)); set(gca, 'YTickLabel',linspace(0,20,11));
set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])
set(gca, 'YDir','reverse'); % reverse y-axis
set(gca, 'Box','on');
xlim([50 1000]);
ylim([60 1950]);
% xlim([1 img_size(2)]);
% ylim([1 img_size(1)]);



%% Make direction movie
% Make figure handle
fh = figure; clf;
% set position on screen
set(fh,'position',[-1850 570 560 420]);
%subplot(1,2,1);

% Show scatterplot with non-zero positions
fig = scatter(0, 0);

% Scatterplot settings
set(fig,'SizeData', 0.5); % size of dots
set(fig,'MarkerFacecolor','flat'); % appearance of dots
xlabel('Lateral [mm]'); ylabel('Axial [mm]'); % title('Micro-Bubble image');
set(gca,'Xtick',linspace(50,1000,6)); set(gca, 'XTickLabel',linspace(0,10,6));
set(gca,'Ytick',linspace(60,1950,11)); set(gca, 'YTickLabel',linspace(0,20,11));
set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])
set(gca, 'YDir','reverse'); % reverse y-axis
set(gca, 'Box','on');
xlim([50 1000]);
ylim([60 1950]);


outputVideo=VideoWriter(['mb_tracking_vid_inflow_1']);
outputVideo.FrameRate=20;
open(outputVideo);
mov(1:size(MB_index_filter_dir,2))= struct('cdata',[],'colormap',[]);
pause(0.2)

i_end = 1;
for i = 1:size(MB_index_filter_dir,2)
    MB_index = MB_index_filter_dir(i);
    
    i_start = i_end;
    i_end = i_end+MB_log_copy(MB_index).age(3)-3;
    
    fig.XData(i_start:i_end) =  MB_log_copy(MB_index).centroid(2:end-1,2)';%scatter(vel_abs_img_list_x(i), vel_abs_img_list_y(i));
    fig.YData(i_start:i_end) =  MB_log_copy(MB_index).centroid(2:end-1,1)';
    fig.CData(i_start:i_end,:) = dir_color(floor(mod(atan2(MB_log_copy(MB_index).vel(2:end,1),MB_log_copy(MB_index).vel(2:end,2))+5/2*pi,2*pi)/(2*pi)*63+1),:);%MB_vel_list_dir(i_start:i_end);% setting colors of individual dots depending on direction
    mov=getframe(gcf);
    writeVideo(outputVideo,mov.cdata);
end

close(gcf)
close(outputVideo);

%% Make vel image
% Make figure handle
fh = figure; clf;

% Colormap
vel_color = autumn;
vel_color = (vel_color);

% Max and min value for colormap
val_range_max = 1.5;
val_range_min = 0;

% set position on screen
set(fh,'position',[-1850 570 560 420]);
% subplot(1,2,1);

% Show scatterplot with non-zero positions
fig = scatter(0, 0);
% fig = scatter(vel_abs_img_list_x, vel_abs_img_list_y);
% fig.CData = dir_color(vel_dir_img_list_val,:); % setting colors of individual dots depending on direction
MB_XData = [];
MB_YData = [];
MB_VelData = [];

i_end = 1;
for i = 1:size(MB_index_filter_dir,2)
    MB_index = MB_index_filter_dir(i);
    
    i_start = i_end;
    i_end = i_end+MB_log_copy(MB_index).age(3)-3;
    
    MB_XData = [MB_XData MB_log_copy(MB_index).centroid(2:end-1,2)'];%scatter(vel_abs_img_list_x(i), vel_abs_img_list_y(i));
    MB_YData = [MB_YData MB_log_copy(MB_index).centroid(2:end-1,1)'];
    MB_VelData = [MB_VelData sqrt((MB_log_copy(MB_index).vel(2:end,2)).^2 + (MB_log_copy(MB_index).vel(2:end,1)).^2)'];
end

MB_VelData = MB_VelData * 10*10^(-6)*fps*10^3;

MB_VelData(find(MB_VelData > val_range_max)) = val_range_max; 
MB_VelData(find(MB_VelData < val_range_min)) = val_range_min;
MB_VelData = MB_VelData-val_range_min;
MB_VelData = floor(MB_VelData/(val_range_max-val_range_min)*(size(vel_color,1)-1)+1);
MB_CData = vel_color(MB_VelData,:);
fig.XData = MB_XData;
fig.YData = MB_YData;
fig.CData = MB_CData;

% Scatterplot settings
set(fig,'SizeData', 0.5); % size of dots
set(fig,'MarkerFacecolor','flat'); % appearance of dots
xlabel('Lateral [mm]'); ylabel('Axial [mm]'); % title('Micro-Bubble image');
set(gca,'Xtick',linspace(50,1000,6)); set(gca, 'XTickLabel',linspace(0,10,6));
set(gca,'Ytick',linspace(60,1950,11)); set(gca, 'YTickLabel',linspace(0,20,11));
set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])
set(gca, 'YDir','reverse'); % reverse y-axis
set(gca, 'Box','on');
xlim([50 1000]);
ylim([60 1950]);
% xlim([1 img_size(2)]);
% ylim([1 img_size(1)]);

% %Colorbar prop
ch = colorbar;
% set(ch,'position',[0.5729 0.5795 0.0306 0.3167]);
 set(ch,'TicksMode','manual');
 set(ch,'Ticks',linspace(0, 1,6));
 set(ch,'TickLabelsMode','manual');
 set(ch,'TickLabels',fliplr(linspace(val_range_min, val_range_max,6)));
 colormap(flipud(vel_color));
 ylabel(ch,'Velocity [mm/s]')


%% Color-wheel
% HSV color circle -- start
% Compute HSV image
rows = 400;
columns = 400;
midX = columns / 2;
midY = rows / 2;
% Construct v image as uniform.
v = 1 * ones(rows, columns);
s = zeros(size(v)); % Initialize.
h = zeros(size(v)); % Initialize.
% Construct the h image as going from 0 to 1 as the angle goes from 0 to 360.
% Construct the S image going from 0 at the center to 1 at the edge.
for c = 1 : columns
	for r = 1 : rows
		% Radius goes from 0 to 1 at edge, a little more in the corners.
		radius = sqrt((r - midY)^2 + (c - midX)^2) / min([midX, midY]);
		s(r, c) = min(1, radius)+0.4; % Max out at 1
		h(r, c) = mod(atan2d((c - midX),(r - midY))+270+color_rotate/(2*pi)*360,360)-180;
    end
end
% Flip h right to left.
h = fliplr(mat2gray(h));
% Construct the hsv image.
hsvImage = cat(3, h, s, v);

% Construct the RGB image.
rgbImage = hsv2rgb(hsvImage);

hold on
% Position Colorwheel
ac = axes('Position',[0.330 0.1320 0.153 0.153]);

% Show colorwheel
imshow(rgbImage);
% Position Colorwheel
ac = axes('Position',[0.770 0.1320 0.153 0.153]);

% Show colorwheel
imshow(rgbImage);
% HSV color circle -- end


%% Density image from scatter_matrix
norm = max(scatter_matrix(:));
log_scatter_matrix = 20*log10(scatter_matrix/norm);
figure(); 
imagesc(log_scatter_matrix,[-40 0]); 
colormap('red');


%% Density image using hist3
hist_data = hist3([subscript_x; subscript_y]',[1000 600]);
hist_data = hist_data';
norm = max(hist_data(:));
log_hist_data = 20*log10(hist_data/norm);
figure();
imagesc(log_hist_data);
colormap(hot);
set(gca, 'YDir','reverse')
set(gca,'Xtick',linspace(0,350,5)); set(gca, 'XTickLabel',linspace(0,12,5));
set(gca,'Ytick',linspace(0,350,6)); set(gca, 'YTickLabel',linspace(0,25,6));

xlabel('Lateral [mm]'); ylabel('Axial [mm]'); % title('Micro-Bubble image');
% set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
% set(gca, 'PlotBoxAspectRatio',[1 1 1])
set(gca, 'Box','on');


%% Paper plot
fh = figure;
% set position on screen
set(fh,'position',[-1850 570 560 420]);
set(fh,'color','w');
% temp = subplot(2,1,1);
set(gca,'Xtick',linspace(235,295,7)); set(gca, 'XTickLabel',linspace(0,600,7));
set(gca,'Ytick',linspace(89,119,7)); set(gca, 'YTickLabel',linspace(0,300,7));
xlabel('Lateral [µm]'); ylabel('Axial [µm]'); % title('Micro-Bubble image');
xlim([235 300]);
ylim([89 118]);
set(gca, 'YDir','reverse'); % reverse y-axis
set(gca, 'Box','on');
set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])

%set(gca, 'XGrid','on')
%set(gca, 'YGrid','on')
%set(gca, 'GridColor',[1 1 1])
%set(gca, 'GridAlpha',0.5)
%set(gca,'Color',[0 0 0])
std_factor = 1;
hold on
%% confidens 1
i = 1;
p1 = polyfit(MB_plot(i).centroid(:,2),MB_plot(i).centroid(:,1),1);
curve_fit1 = polyval(p1,MB_plot(i).centroid(:,2));

var_curve1 = sum((curve_fit1-MB_plot(i).centroid(:,1)).^2/(size(curve_fit1,1)-1));
std_curve1 = sqrt(var_curve1);

p1 = polyfit(MB_plot(i).centroid(1:120,2),MB_plot(i).centroid(1:120,1),1);

curve_x1 = min(MB_plot(i).centroid(1:120,2)):max(MB_plot(i).centroid(1:120,2));
curve_mid1 = polyval(p1,curve_x1);
%curve_mid_vector = [min(MB_plot(i).centroid(:,2)):max(MB_plot(i).centroid(:,2)); curve_mid]; 
%curve_top_vector = curve_mid_vector(:,1:end-1) + (curve_mid_vector(:,2:end)-curve_mid_vector(:,1:end-1))./[sqrt((curve_mid_vector(1,2:end)-curve_mid_vector(1,1:end-1)).^2+(curve_mid_vector(2,2:end)-curve_mid_vector(2,1:end-1)).^2); sqrt((curve_mid_vector(1,2:end)-curve_mid_vector(1,1:end-1)).^2+(curve_mid_vector(2,2:end)-curve_mid_vector(2,1:end-1)).^2)].*[ones(1,size(curve_mid_vector,2)-1); ones(1,size(curve_mid_vector,2)-1)]*(3*std_curve);
curve_top1 = curve_mid1 + std_factor * std_curve1;
curve_bot1 = curve_mid1 - std_factor * std_curve1;

curve_fill_x1 = [curve_x1, fliplr(curve_x1)];
curve_fill_y1 = [curve_top1, fliplr(curve_bot1)];

% figure();
fill_h1 = fill(curve_fill_x1,curve_fill_y1,'g');
% hold on
% plot(curve_x1,curve_mid1,'r-');
% plot(curve_x1,curve_top1,'g-');
% plot(curve_x1,curve_mid1,'g--');
% plot(curve_x1,curve_bot1,'g-');

%plot( MB_plot(i).centroid(:,2), MB_plot(i).centroid(:,1), 'bo');
%ax_h = gca;
set(fill_h1,'facealpha',.1)
set(fill_h1,'EdgeColor','None');
%set(ax_h,'YDir','reverse');

%% confidens 2
i = 3;
p2 = polyfit(MB_plot(i).centroid(:,2),MB_plot(i).centroid(:,1),1);
curve_fit2 = polyval(p2,MB_plot(i).centroid(:,2));

var_curve2 = sum((curve_fit2-MB_plot(i).centroid(:,1)).^2/(size(curve_fit2,1)-1));
std_curve2 = sqrt(var_curve2);

p2 = polyfit(MB_plot(i).centroid(5:42,2),MB_plot(i).centroid(5:42,1),1);

curve_x2 = min(MB_plot(i).centroid(5:42,2)):max(MB_plot(i).centroid(5:42,2));
curve_mid2 = polyval(p2,curve_x2);
%curve_mid_vector = [min(MB_plot(i).centroid(:,2)):max(MB_plot(i).centroid(:,2)); curve_mid]; 
%curve_top_vector = curve_mid_vector(:,1:end-1) + (curve_mid_vector(:,2:end)-curve_mid_vector(:,1:end-1))./[sqrt((curve_mid_vector(1,2:end)-curve_mid_vector(1,1:end-1)).^2+(curve_mid_vector(2,2:end)-curve_mid_vector(2,1:end-1)).^2); sqrt((curve_mid_vector(1,2:end)-curve_mid_vector(1,1:end-1)).^2+(curve_mid_vector(2,2:end)-curve_mid_vector(2,1:end-1)).^2)].*[ones(1,size(curve_mid_vector,2)-1); ones(1,size(curve_mid_vector,2)-1)]*(3*std_curve);
curve_top2 = curve_mid2 + std_factor * std_curve2;
curve_bot2 = curve_mid2 - std_factor * std_curve2;

curve_fill_x2 = [curve_x2, fliplr(curve_x2)];
curve_fill_y2 = [curve_top2, fliplr(curve_bot2)];

% figure();
fill_h2 = fill(curve_fill_x2,curve_fill_y2, 'b');
% hold on
%plot(curve_x,curve_mid,'r-');
% plot(curve_x2,curve_top2,'b-');
% plot(curve_x2,curve_mid2,'b--');
% plot(curve_x2,curve_bot2,'b-');

% plot( MB_plot(i).centroid(:,2), MB_plot(i).centroid(:,1), 'bo');
%ax_h = gca;
set(fill_h2,'facealpha',.1)
set(fill_h2,'EdgeColor','None');
%set(ax_h,'YDir','reverse');


%%
% Show scatterplot with non-zero positions

fig = scatter(vel_abs_img_list_x, vel_abs_img_list_y);



% Scatterplot settings
set(fig,'SizeData', 10); % size of dots
set(fig,'MarkerFacecolor','flat'); % appearance of dots
fig.CData = dir_color(vel_dir_img_list_val,:); % setting colors of individual dots depending on direction

plot(curve_x1,curve_top1,'g-');
plot(curve_x1,curve_mid1,'g--');
plot(curve_x1,curve_bot1,'g-');

plot(curve_x2,curve_top2,'b-');
plot(curve_x2,curve_mid2,'b--');
plot(curve_x2,curve_bot2,'b-');

%% video

fig = scatter(0, 0);
% Scatterplot settings
set(fig,'SizeData', 10); % size of dots
set(fig,'MarkerFacecolor','flat'); % appearance of dots

MB_XData = [];
MB_YData = [];
MB_CData = [];

i_end = 1;
for i = 1:size(MB_index_filter_dir,2)
    MB_index = MB_index_filter_dir(i);
    
    i_start = i_end;
    i_end = i_end+MB_log_copy(MB_index).age(3)-3;
    
    MB_XData(i_start:i_end) =  MB_log_copy(MB_index).centroid(2:end-1,2)';%scatter(vel_abs_img_list_x(i), vel_abs_img_list_y(i));
    MB_YData(i_start:i_end) =  MB_log_copy(MB_index).centroid(2:end-1,1)';
    MB_CData(i_start:i_end,:) = dir_color(floor(mod(atan2(MB_log_copy(MB_index).vel(2:end,1),MB_log_copy(MB_index).vel(2:end,2))+5/2*pi,2*pi)/(2*pi)*63+1),:);%MB_vel_list_dir(i_start:i_end);% setting colors of individual dots depending on direction 
end

outputVideo=VideoWriter(['mb_tracking_vid_mc']);
outputVideo.FrameRate=15;
open(outputVideo);
mov(1:size(MB_index_filter_dir,2))= struct('cdata',[],'colormap',[]);
pause(0.2)

for i = 1:size(MB_XData,2)
    
    fig.XData(i) =  MB_XData(i);
    fig.YData(i) =  MB_YData(i);
    fig.CData(i,:) = MB_CData(i,:)
    mov=getframe(gcf);
    writeVideo(outputVideo,mov.cdata);
end

close(gcf)
close(outputVideo);

%% Colorwheel
% Position Colorwheel
ac = axes('Position',[0.8 0.08 0.153 0.153]);

% Show colorwheel
imshow(rgbImage);
% Position Colorwheel
ac = axes('Position',[0.8 0.545 0.153 0.153]);

% Show colorwheel
imshow(rgbImage);
% HSV color circle -- end



%% B-mode / Doppler image
load('/data/cfudata6/s134082/micro_bubble_data/mat_files/console_video/video_console.mat')
figure();
subplot(1,2,1);
mov_b_frame = mov_b(:,:,:,30);
imagesc(mov_b_frame(20:600,100:380,:));
rectangle('Position',[50 300 30 10],'EdgeColor','r','LineWidth',2);
annotation('textarrow',[0.31 0.24],[0.6 0.51],'color','r','LineWidth',2)
set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])
set(gca,'Xtick',linspace(1,280,6)); set(gca, 'XTickLabel',linspace(0,10,6));
set(gca,'Ytick',linspace(1,580,11)); set(gca, 'YTickLabel',linspace(0,20,11));
xlabel('Lateral [mm]'); ylabel('Axial [mm]'); % title('Micro-Bubble image');

subplot(1,2,2)
mov_dob_frame = mov_dob(:,:,:,12);
mov_dob_frame(23:144,360:374,:) = mov_dob_frame(48:169,863:877,:);

imagesc(mov_dob_frame(20:600,100:380,:));
set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])
set(gca,'Xtick',linspace(1,280,6)); set(gca, 'XTickLabel',linspace(0,10,6));
set(gca,'Ytick',linspace(1,580,11)); set(gca, 'YTickLabel',linspace(0,20,11));
xlabel('Lateral [mm]'); ylabel('Axial [mm]'); % title('Micro-Bubble image');



%% Motion tracking
figure();
% subplot(1,2,1);
load('/data/cfudata6/s134082/MB_tracking_matlab/Movement_compensation/MT_with_pulse_3_rep.mat')
plot(axial_vel_mean(44,:),'k');
xlim([1 1290]); ylim([-2 2.5]);
set(gca,'Xtick',linspace(1,1290,8)); set(gca, 'XTickLabel',linspace(0,42,8));
set(gca,'Ytick',linspace(-2,2,5)); set(gca, 'YTickLabel',linspace(-50,50,5));
xlabel('Time [ms]'); ylabel('Axial [µm]'); % title('Micro-Bubble image');

%%
% subplot(1,2,2)
figure();
load('/data/cfudata6/s134082/MB_tracking_matlab/Movement_compensation/MT_without_pulse.mat')
plot(axial_vel_mean(44,:),'k');
xlim([1 430]); ylim([-2 2.5]);
set(gca,'Xtick',linspace(1,370,3)); set(gca, 'XTickLabel',linspace(0,12,3));
set(gca,'Ytick',linspace(-2,2,5)); set(gca, 'YTickLabel',linspace(-50,50,5));
xlabel('Time [ms]'); ylabel('Axial [µm]'); % title('Micro-Bubble image');

%%
figure();
load('/data/cfudata6/s134082/MB_tracking_matlab/Movement_compensation/MT_pulse_1_rep.mat')
plot(axial_vel_mean(44,:),'k');
xlim([1 430]); ylim([-2 2.5]);
set(gca,'Xtick',linspace(1,370,3)); set(gca, 'XTickLabel',linspace(0,12,3));
set(gca,'Ytick',linspace(-2,2,5)); set(gca, 'YTickLabel',linspace(-50,50,5));
xlabel('Time [ms]'); ylabel('Axial [µm]'); % title('Micro-Bubble image');
