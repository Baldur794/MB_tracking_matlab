%% Axes for plots

% Get img size
img_size = size(load_img(1, img_type, folderName, img_resolution, area_of_interest)); % [y,x]

% Propagation speed of sound
c = 1540;

% Create axes
t_axis = linspace(area_of_interest.axial_init*(1/sampling_rate),area_of_interest.axial_end*(1/sampling_rate),img_size(1));
depth_axis = (t_axis*c)/2;   %  Values always in MKS :end_line)
lateral_axis = linspace((area_of_interest.lateral_init-1)*element_pitch/line_density,(area_of_interest.lateral_end-1)*element_pitch/line_density,img_size(2));
% lateral_axis = linspace((area_of_interest.lateral_init-floor((area_of_interest.lateral_end-(area_of_interest.lateral_init-1))/2))*element_pitch/line_density,(area_of_interest.lateral_end-floor((area_of_interest.lateral_end-(area_of_interest.lateral_init-1))/2))*element_pitch/line_density,img_size_lateral);

% Frame rate
fps = 163;


%%
MB_log_copy = MB_log; % Copy of MB_log
MB_index_list = 1:size(MB_log_copy,2); % List of indexes for valid MBs

% Setup conditions
% Use only MBs found within min and max frame number condition
MB_frameNumber_min_condition_flag = false;
MB_frameNumber_min_condition = 550;
MB_frameNumber_max_condition = 590;
MB_index_min_max_frameNumber = [];

% Use only MBs with an age between min and max age condition
MB_age_min_max_condition_flag = true;
MB_age_min_condition = 6;
MB_age_max_condition = 3000;
MB_index_min_max_age = [];

% Use only MBs with a blob count below max condition
MB_max_blob_count_condition_flag = false;
MB_max_blob_count_condition = 20; 
MB_index_max_blob_count = [];

% Use only MBs with an avg velocity above min condition
MB_min_avg_vel_condition_flag = true;
MB_min_avg_vel_condition = 1;
MB_index_min_avg_vel = []; 

% Use only MBs located around other MBs with conditions about distance and avg density.
MB_min_avg_density_condition_flag = false;
MB_min_avg_density_condition = 0.5; 
MB_window_size_density_avg = [5 5]; % Window for density condition
MB_index_min_avg_density = []; 

% Use only MBs for determine inflow/outflow. 0 = not used
inflow_outflow_center_coord_condition_flag = false;
inflow_outflow_const = 1; % Towards center coord set 1. Away set -1.
inflow_outflow_center_coord_condition = [0 0]; % [y,x]
MB_index_inflow_outflow = []; % Filtered list of MBs satisfying conditions

% Smooth velocity/direction for each MB
smooth_vel_dir_flag = true;
smooth_vel_dir_length = 50; % Number of MBs used for smoothing

%
% Check frame condition
if MB_frameNumber_min_condition_flag; 
    for i = 1:size(MB_index_list,2)
        MB_index = MB_index_list(i);
        % Check condition
        if (MB_log_copy(MB_index).age(1) < MB_frameNumber_min_condition) || (MB_log_copy(MB_index).age(2) > MB_frameNumber_max_condition)
            MB_index_list(i) = 0;
        end
    end
    MB_index_list(MB_index_list == 0) = [];
    MB_index_min_max_frameNumber = MB_index_list;
end


% Check age condition
if MB_age_min_max_condition_flag;
    for i = 1:size(MB_index_list,2)
        MB_index = MB_index_list(i);
        % Check condition
        if (MB_log_copy(MB_index).age(3) < MB_age_min_condition) || (MB_log_copy(MB_index).age(3) > MB_age_max_condition)
            MB_index_list(i) = 0;
        end
    end
    MB_index_list(MB_index_list == 0) = [];
    MB_index_min_max_age = MB_index_list;
end


% Check max blob count condition
if MB_max_blob_count_condition_flag;
    for i = 1:size(MB_index_list,2)
        MB_index = MB_index_list(i);
        % Check condition
        if (max(MB_log_copy(MB_index).count(:)) > MB_max_blob_count_condition)
            MB_index_list(i) = 0;
        end
    end
    MB_index_list(MB_index_list == 0) = [];
    MB_index_max_blob_count = MB_index_list;
end


% Check min average speed condition
if MB_min_avg_vel_condition_flag
    for i = 1:size(MB_index_list,2)
        MB_index = MB_index_list(i);
        avg_vel = sqrt((MB_log_copy(MB_index).centroid(1,1)-MB_log_copy(MB_index).centroid(end,1))^2+(MB_log_copy(MB_index).centroid(1,2)-MB_log_copy(MB_index).centroid(end,2))^2)/MB_log_copy(MB_index).age(3);      
        % Check condition
        if (avg_vel < MB_min_avg_vel_condition)
            MB_index_list = 0;
        end
    end
    MB_index_list(MB_index_list == 0) = [];
    MB_index_min_avg_vel = MB_index_list;
end

% Check avg density condition
if MB_min_avg_density_condition_flag
    % Scatter image from all MB locations
    scatter_matrix = zeros(img_size);
    for i = 1:size(MB_index_list ,2)
        MB_index = MB_index_list(i);
        scatter_matrix(sub2ind(img_size,MB_log_copy(MB_index).centroid(1:end-1,1) ,MB_log_copy(MB_index).centroid(1:end-1,2))) = scatter_matrix(sub2ind(img_size,MB_log_copy(MB_index).centroid(1:end-1,1) ,MB_log_copy(MB_index).centroid(1:end-1,2))) +1;
    end
    
    for i = 1:size(MB_index_list ,2)
        MB_index = MB_index_list(i);
        
        MB_density_list = []; % Contains number of neighbouring MBs around each MB position
        MB_length = size(MB_log_copy(MB_index).centroid(1:end-1,1),1);
        
        % Calculated density of MB in window
        for j = 1:MB_length
            % Update window coordinates
            MB_window_coord(1,1) = MB_log_copy(MB_index).centroid(j,1)-MB_window_size_density_avg(1); % y1
            MB_window_coord(2,1) = MB_log_copy(MB_index).centroid(j,1)+MB_window_size_density_avg(1); % y2
            MB_window_coord(1,2) = MB_log_copy(MB_index).centroid(j,2)-MB_window_size_density_avg(2); % x1
            MB_window_coord(2,2) = MB_log_copy(MB_index).centroid(j,2)+MB_window_size_density_avg(2); % x2
            
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
            MB_density_list(j) = sum(sum(scatter_matrix(MB_window_coord(1,1):MB_window_coord(2,1),MB_window_coord(1,2):MB_window_coord(2,2))))-1; % number of MBs found around given MB with # MB_index
        end
        % Calculate average
        MB_density_sum = sum(MB_density_list)/MB_length;
        
        % Check density condition
        if (MB_density_sum < MB_min_avg_density_condition)% && (min(MB_density_list) >= MB_index_inflow_outflow)
            MB_index_list(i) = 0;
        end
    end
    MB_index_list(MB_index_list == 0) = [];
    MB_index_min_avg_density = MB_index_list;
end

% Find flow towards/away from center coord
if inflow_outflow_center_coord_condition_flag
    for i = 1:size(MB_index_list ,2)
        MB_index = MB_index_list(i);   
        distance_to_center_start = sqrt((MB_log_copy(MB_index).centroid(1,1)-inflow_outflow_center_coord_condition(1))^2+(MB_log_copy(MB_index).centroid(1,2)-inflow_outflow_center_coord_condition(2))^2);
        distance_to_center_end = sqrt((MB_log_copy(MB_index).centroid(end,1)-inflow_outflow_center_coord_condition(1))^2+(MB_log_copy(MB_index).centroid(end,2)-inflow_outflow_center_coord_condition(2))^2);
        % Check density condition
        if ((distance_to_center_start - distance_to_center_end)*inflow_outflow_const)
            MB_index_list = 0;
        end
    end
    MB_index_list(MB_index_list == 0) = [];
    MB_index_inflow_outflow = MB_index_list;
end

% Smooth velocities
if smooth_vel_dir_flag
    if smooth_vel_dir_length > MB_age_min_condition-2
        smooth_vel_dir_length = MB_age_min_condition-2;
    end
    for i = 1:size(MB_index_list,2)
        MB_index = MB_index_list(i);
        for j = 1:MB_log_copy(MB_index).age(3)-1
            if floor(j-(smooth_vel_dir_length-1)/2) < 1
                MB_log_copy(MB_index).vel(j,1) = sum(MB_log_copy(MB_index).vel(1:1+(smooth_vel_dir_length-1),1))/smooth_vel_dir_length;
                MB_log_copy(MB_index).vel(j,2) = sum(MB_log_copy(MB_index).vel(1:1+(smooth_vel_dir_length-1),2))/smooth_vel_dir_length;
            elseif floor(j+(smooth_vel_dir_length-1)/2) > MB_log_copy(MB_index).age(3)-1
                MB_log_copy(MB_index).vel(j,1) = sum(MB_log_copy(MB_index).vel(end-(smooth_vel_dir_length-1):end-1,1))/smooth_vel_dir_length;
                MB_log_copy(MB_index).vel(j,2) = sum(MB_log_copy(MB_index).vel(end-(smooth_vel_dir_length-1):end,2))/smooth_vel_dir_length;
            else
                MB_log_copy(MB_index).vel(j,1) = sum(MB_log_copy(MB_index).vel(floor(j-(smooth_vel_dir_length-1)/2):floor(j+(smooth_vel_dir_length-1)/2),1))/smooth_vel_dir_length;
                MB_log_copy(MB_index).vel(j,2) = sum(MB_log_copy(MB_index).vel(floor(j-(smooth_vel_dir_length-1)/2):floor(j+(smooth_vel_dir_length-1)/2),2))/smooth_vel_dir_length;
            end
        end
    end
end
 

% Make direction image
% Make figure handle
fh = figure; clf;
% set position on screen
set(fh,'position',[-1850 570 560 420]);
% subplot(1,2,1);

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
for i = 1:size(MB_index_list,2)
    MB_index = MB_index_list(i);
    
    i_start = i_end;
    i_end = i_end+MB_log_copy(MB_index).age(3)-3;
    
    MB_XData = [MB_XData MB_log_copy(MB_index).centroid(2:end-1,2)'];%scatter(vel_abs_img_list_x(i), vel_abs_img_list_y(i));
    MB_YData = [MB_YData MB_log_copy(MB_index).centroid(2:end-1,1)'];
    MB_CData = [MB_CData; dir_color(floor(mod(atan2(MB_log_copy(MB_index).vel(2:end,1),MB_log_copy(MB_index).vel(2:end,2))+5/2*pi,2*pi)/(2*pi)*63+1),:)];% setting colors of individual dots depending on direction
end
% MB_index_list = [];
% MB_coord_list = [];
% MB_vel_list = [];
fig.XData = MB_XData;
fig.YData = MB_YData;
fig.CData = MB_CData;

% Scatterplot settings
set(fig,'SizeData', 20);%0.5); % size of dots
set(fig,'MarkerFacecolor','flat'); % appearance of dots
xlabel('Lateral [mm]'); ylabel('Axial [mm]'); % title('Micro-Bubble image');
set(gca,'Xtick',linspace(0,size(lateral_axis,2),6)); set(gca, 'XTickLabel',linspace(round(lateral_axis(1)*1000),round(lateral_axis(end)*1000),6));
set(gca,'Ytick',linspace(0,size(depth_axis,2),6)); set(gca, 'YTickLabel',linspace(round(depth_axis(1)*1000),round(depth_axis(end)*1000),6));
set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])
set(gca, 'YDir','reverse'); % reverse y-axis
set(gca, 'Box','on');
xlim([0 size(lateral_axis,2)]);
ylim([0 size(depth_axis,2)]);


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
set(gca,'Xtick',linspace(0,945,6)); set(gca, 'XTickLabel',linspace(lateral_axis(1),lateral_axis(end),6));
set(gca,'Ytick',linspace(0,769,6)); set(gca, 'YTickLabel',linspace(depth_axis(1),depth_axis(end),6));
set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])
set(gca, 'YDir','reverse'); % reverse y-axis
set(gca, 'Box','on');
% xlim([50 1000]);
% ylim([60 1950]);


outputVideo=VideoWriter(['mb_tracking_vid_inflow_1']);
outputVideo.FrameRate=20;
open(outputVideo);
mov(1:size(MB_index_list,2))= struct('cdata',[],'colormap',[]);
pause(0.2)

i_end = 1;
for i = 1:size(MB_index_list,2)
    MB_index = MB_index_list(i);
    
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
for i = 1:size(MB_index_list,2)
    MB_index = MB_index_list(i);
    
    i_start = i_end;
    i_end = i_end+MB_log_copy(MB_index).age(3)-3;
    
    MB_XData = [MB_XData MB_log_copy(MB_index).centroid(2:end-1,2)'];%scatter(vel_abs_img_list_x(i), vel_abs_img_list_y(i));
    MB_YData = [MB_YData MB_log_copy(MB_index).centroid(2:end-1,1)'];
    MB_VelData = [MB_VelData sqrt((MB_log_copy(MB_index).vel(2:end,2)).^2 + (MB_log_copy(MB_index).vel(2:end,1)).^2)'];
end

% Max and min value for colormap
val_range_max = 25;
val_range_min = 5;

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
set(fig,'SizeData', 20);%0.5); % size of dots
set(fig,'MarkerFacecolor','flat'); % appearance of dots
xlabel('Lateral [mm]'); ylabel('Axial [mm]'); % title('Micro-Bubble image');
set(gca,'Xtick',linspace(0,945,6)); set(gca, 'XTickLabel',linspace(round(lateral_axis(1)*1000),round(lateral_axis(end)*1000),6));
set(gca,'Ytick',linspace(0,769,6)); set(gca, 'YTickLabel',linspace(round(depth_axis(1)*1000),round(depth_axis(end)*1000),6));
set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])
set(gca, 'YDir','reverse'); % reverse y-axis
set(gca, 'Box','on');

%Colorbar prop
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
% % Position Colorwheel
% ac = axes('Position',[0.330 0.1320 0.153 0.153]);
% 
% % Show colorwheel
% imshow(rgbImage);
% Position Colorwheel
ac = axes('Position',[0.350 0.24 0.153 0.153]);

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
hist_data = hist3(fliplr(MB_coord_list'),[1000 600]);
hist_data = hist_data';
norm = max(hist_data(:));
log_hist_data = 20*log10(hist_data/norm);
figure();
imagesc(log_hist_data);
colormap(hot);
% set(gca, 'YDir','reverse')
set(gca,'Xtick',linspace(0,350,5)); set(gca, 'XTickLabel','');
set(gca,'Ytick',linspace(0,350,6)); set(gca, 'YTickLabel','');

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
for i = 1:size(MB_index_list,2)
    MB_index = MB_index_list(i);
    
    i_start = i_end;
    i_end = i_end+MB_log_copy(MB_index).age(3)-3;
    
    MB_XData(i_start:i_end) =  MB_log_copy(MB_index).centroid(2:end-1,2)';%scatter(vel_abs_img_list_x(i), vel_abs_img_list_y(i));
    MB_YData(i_start:i_end) =  MB_log_copy(MB_index).centroid(2:end-1,1)';
    MB_CData(i_start:i_end,:) = dir_color(floor(mod(atan2(MB_log_copy(MB_index).vel(2:end,1),MB_log_copy(MB_index).vel(2:end,2))+5/2*pi,2*pi)/(2*pi)*63+1),:);%MB_vel_list_dir(i_start:i_end);% setting colors of individual dots depending on direction 
end

outputVideo=VideoWriter(['mb_tracking_vid_mc']);
outputVideo.FrameRate=15;
open(outputVideo);
mov(1:size(MB_index_list,2))= struct('cdata',[],'colormap',[]);
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
set(gca,'Xtick',linspace(1,1290,6)); set(gca, 'XTickLabel',linspace(0,2.5,6));
set(gca,'Ytick',linspace(-2,2,5)); set(gca, 'YTickLabel',linspace(-50,50,5));
xlabel('Time [s]'); ylabel('Axial [µm]'); % title('Micro-Bubble image');

%%
% subplot(1,2,2)
figure();
load('/data/cfudata6/s134082/MB_tracking_matlab/Movement_compensation/MT_without_pulse.mat')
plot(axial_vel_mean(44,:),'k');
xlim([1 430]); ylim([-2 2.5]);
set(gca,'Xtick',linspace(1,370,3)); set(gca, 'XTickLabel',linspace(0,0.8,3));
set(gca,'Ytick',linspace(-2,2,5)); set(gca, 'YTickLabel',linspace(-50,50,5));
xlabel('Time [s]'); ylabel('Axial [µm]'); % title('Micro-Bubble image');

%%
figure();
load('/data/cfudata6/s134082/MB_tracking_matlab/Movement_compensation/MT_pulse_1_rep.mat')
plot(axial_vel_mean(44,:),'k');
xlim([1 430]); ylim([-2 2.5]);
set(gca,'Xtick',linspace(1,370,3)); set(gca, 'XTickLabel',linspace(0,0.8,3));
set(gca,'Ytick',linspace(-2,2,5)); set(gca, 'YTickLabel',linspace(-50,50,5));
xlabel('Time [ms]'); ylabel('Axial [µm]'); % title('Micro-Bubble image');
