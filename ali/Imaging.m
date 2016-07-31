%% Age and Density filter
MB_index_list = []; % List of indexes with for found MB's
MB_vel_list = []; % List with all velocity vectors
MB_log_copy = MB_log;
% List all MB positions and their velocities

for MB_index = 1:size(MB_log_copy,2)
   if (MB_log_copy(MB_index).age(3) > MB_age_condition) && (max(MB_log_copy(MB_index).count(:)) <= MB_count_condition)
       MB_index_list = [MB_index_list, sub2ind([img_size(1:2)],MB_log_copy(MB_index).centroid(1:end-1,1) ,MB_log_copy(MB_index).centroid(1:end-1,2))'];
       MB_vel_list = [MB_vel_list, [MB_log_copy(MB_index).centroid(2:end,1)-MB_log_copy(MB_index).centroid(1:end-1,1),MB_log_copy(MB_index).centroid(2:end,2)-MB_log_copy(MB_index).centroid(1:end-1,2)]'];
   end
end

% Scatter image from all MB locations
scatter_matrix = zeros([img_size(1:2)]); 
for i = 1:size(MB_index_list,2)
    scatter_matrix(MB_index_list(i)) = scatter_matrix(MB_index_list(i)) + 1;
end

% Velocity img, Amplitude and direction
% Final velocity and direction img
vel_abs_img = zeros([img_size(1:2)]);
vel_dir_img = zeros([img_size(1:2)]);

% Weighing filter for weighted velocities
%weighing_filter_radius = 3; % Radius around centroid to be considered
weighing_filter_size = 2*weighing_filter_radius+1; % Must be an odd number
weighing_filter_halfsize = floor(weighing_filter_size/2); % Only to make indexing easier
weighing_filter = zeros(weighing_filter_size); % Contains the final filter
weighing_center = [ceil(weighing_filter_size/2),ceil(weighing_filter_size/2)]; % Center coordinates
%weighing_factor = 1; % Distance weighing factor

% Makes linear weighing filter for full square
for i = 1:size(weighing_filter,1)
    for j = 1:size(weighing_filter,2)
        weighing_filter(i,j) = sqrt((i-weighing_center(1))^2+(j-weighing_center(2))^2);
    end
end

% Discard indexes with larger radius than weighing_filter_radius
if weighing_factor == 0
    weighing_filter = zeros(size(weighing_filter));
    weighing_filter(ceil(size(weighing_filter,1)/2),ceil(size(weighing_filter,2)/2)) = 1;
else
    temp_weighing_filter = weighing_filter;
    temp_weighing_filter(weighing_filter > weighing_filter_radius) = 0;
    temp_weighing_filter(weighing_filter <= weighing_filter_radius) = 1;
    weighing_filter = weighing_filter * 1/weighing_factor;
    weighing_filter = exp(-(weighing_filter/weighing_filter_radius).^2);
    weighing_filter = weighing_filter.*temp_weighing_filter;
end

% Velocity matrix
vel_logical_matrix = zeros([img_size(1:2)]);
vel_x = zeros([img_size(1:2)]); % x-part of velocity
vel_y = zeros([img_size(1:2)]); % y-part of velocity

for i = 1:size(MB_index_list,2)
    vel_x(MB_index_list(i)) = vel_x(MB_index_list(i)) + MB_vel_list(2,i);
    vel_y(MB_index_list(i)) = vel_y(MB_index_list(i)) + MB_vel_list(1,i);
    vel_logical_matrix(MB_index_list(i)) = vel_logical_matrix(MB_index_list(i)) + 1;
end

% Weight multiple occuring positions
vel_x(find(vel_x > 0)) = vel_x(find(vel_x > 0))./vel_logical_matrix(find(vel_x > 0));
vel_y(find(vel_y > 0)) = vel_y(find(vel_y > 0))./vel_logical_matrix(find(vel_y > 0));
vel_logical_matrix(vel_logical_matrix > 1) = 1;

% y,x subscripts from direct matrix indexes
[subscript_y, subscript_x] = ind2sub([img_size(1:2)],MB_index_list);

% Temporary window and coordinates area in interst
vel_window_coord = [];
vel_window = []; 
vel_window_y = [];
vel_window_x = [];

for i = 1:size(MB_index_list,2)
    
    % Update window coordinates
    vel_window_coord(1,1) = subscript_y(i)-weighing_filter_halfsize;
    vel_window_coord(2,1) = subscript_y(i)+weighing_filter_halfsize;
    vel_window_coord(1,2) = subscript_x(i)-weighing_filter_halfsize;
    vel_window_coord(2,2) = subscript_x(i)+weighing_filter_halfsize;
    
    % Check for out of bounce
    if (vel_window_coord(1,1) > 0) && (vel_window_coord(2,1) < img_size(1)) && (vel_window_coord(1,2) > 0) && (vel_window_coord(2,2) < img_size(2))
        % Windows of interest around specific position
        vel_window = vel_logical_matrix(vel_window_coord(1,1):vel_window_coord(2,1),vel_window_coord(1,2):vel_window_coord(2,2));
        vel_window_y = vel_y(vel_window_coord(1,1):vel_window_coord(2,1),vel_window_coord(1,2):vel_window_coord(2,2));
        vel_window_x = vel_x(vel_window_coord(1,1):vel_window_coord(2,1),vel_window_coord(1,2):vel_window_coord(2,2));
        
        % Distance weighting exp(-(d_i/r)^2)
        vel_weight_matrix = weighing_filter.*vel_window;
        
        % Sum of all weightings
        vel_weight_Z = sum(vel_weight_matrix(:)); 
        
        % Final weighted velocities
        vel_weight_y = sum(sum(1/vel_weight_Z*vel_weight_matrix.*vel_window_y));
        vel_weight_x = sum(sum(1/vel_weight_Z*vel_weight_matrix.*vel_window_x));
        
        % insert weighted velocities to right positions
        vel_y(MB_index_list(i)) = vel_weight_y;
        vel_x(MB_index_list(i)) = vel_weight_x;
    end
end

% Calculate absolute velocity
vel_abs_img = sqrt(vel_y.^2+vel_x.^2);
% Correct for pixel resolution 
vel_abs_img = vel_abs_img; %* 10*10^(-6)*fps*10^3;

% Color map for direction img
colormap_vel_dir = hsv;
colormap_vel_dir(1,:) = [0 0 0];

% Direction img. Outputs in radians
vel_dir_img = atan2(vel_y,vel_x);
% From [-pi,pi] -> [0,2*pi]
vel_dir_img(vel_dir_img ~= 0) = vel_dir_img(vel_dir_img ~= 0) + pi;
% Turn colordisk
color_rotate = 3/2*pi;
vel_dir_img(vel_dir_img ~= 0) = mod(vel_dir_img(vel_dir_img ~= 0) + color_rotate,2*pi);

% Subscripts and index of non zero data pixels (both for vel and dir)
[vel_abs_img_list_y vel_abs_img_list_x] = ind2sub(size(vel_abs_img),find(vel_abs_img(:) > 0)); 
vel_abs_img_list = find(vel_abs_img(:));


% Dir scatter plot
% colormapping for direction img
dir_color = hsv; 

% Make figure handle
fh = figure; clf;
% set position on screen
set(fh,'position',[-1850 570 560 420]);

% List dir values at non-zero indexes
vel_dir_img_list_val = vel_dir_img(vel_abs_img_list);
% Adjust dir values to lay between [1:64] for colormap
vel_dir_img_list_val = floor(vel_dir_img_list_val/(2*pi)*63+1);

% Show scatterplot with non-zero positions
fig = scatter(vel_abs_img_list_x, vel_abs_img_list_y);

% Scatterplot settings
set(fig,'SizeData', 3); % size of dots
set(fig,'MarkerFacecolor','flat'); % appearance of dots
fig.CData = dir_color(vel_dir_img_list_val,:); % setting colors of individual dots depending on direction

xlabel('Lateral [mm]'); ylabel('Axial [mm]'); % title('Micro-Bubble image');
% set(gca,'Xtick',linspace(0,1189,5)); set(gca, 'XTickLabel',linspace(0,12,5));
% set(gca,'Ytick',linspace(0,2489,6)); set(gca, 'YTickLabel',linspace(0,25,6));
% set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
% set(gca, 'PlotBoxAspectRatio',[1 1 1])
xlim([1 img_size(2)]);
ylim([1 img_size(1)]);

set(gca, 'YDir','reverse'); % reverse y-axis
set(gca, 'Box','on');



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
ac = axes('Position',[0.4740 0.1120 0.703 0.153]);

% Show colorwheel
imshow(rgbImage);
% HSV color circle -- end
%% Density image
norm = max(scatter_matrix(:));
log_scatter_matrix = 20*log10(scatter_matrix/norm);
figure(); 
imagesc(log_scatter_matrix,[-40 0]); 
colormap('gray');

%%
figure();
hist3([subscript_x; subscript_y]',[130 50])
hold on

hist_data = hist3([subscript_x; subscript_y]',[130 50]);
hist_data(size(hist_data,1) + 1, size(hist_data,2) + 1) = 0;
xb = linspace(min(subscript_x),max(subscript_x),size(hist_data,2));
yb = linspace(min(subscript_y),max(subscript_y),size(hist_data,1));

img_handle = pcolor(xb,yb,hist_data);
img_handle.ZData = ones(size(hist_data))*2*-max(hist_data(:));
colormap(hot);
set(gca, 'YDir','reverse')
view(3);


%%
%hist_data = hist3([subscript_x; subscript_y]',round(size(img)/10));
hist_data = hist3([subscript_x; subscript_y]',[350 350]);
hist_data = hist_data';
norm = max(hist_data(:));
log_hist_data = 20*log10(hist_data/norm);
figure();
imagesc(hist_data, [0 160]);
colormap(hot);
set(gca, 'YDir','reverse')
xlabel('Lateral [mm]'); ylabel('Axial [mm]'); % title('Micro-Bubble image');
set(gca,'Xtick',linspace(0,350,5)); set(gca, 'XTickLabel',linspace(0,12,5));
set(gca,'Ytick',linspace(0,350,6)); set(gca, 'YTickLabel',linspace(0,25,6));

%%
%% Vel scatter Zoom plot
% Colormap
vel_color = autumn;
vel_color = flipud(vel_color);

% Make figure handle
fh = figure(); clf;
% set position on screen
set(fh,'position',[-1850 570 560 420]);

% Max and min value for colormap
val_range_max = 3;
val_range_min = 0;

% Fitting data within max and min into colormap size
vel_abs_img_list_val = vel_abs_img(vel_abs_img_list);
vel_abs_img_list_val(find(vel_abs_img_list_val > val_range_max)) = val_range_max; 
vel_abs_img_list_val(find(vel_abs_img_list_val < val_range_min)) = val_range_min;
vel_abs_img_list_val = vel_abs_img_list_val-val_range_min;
vel_abs_img_list_val = floor(vel_abs_img_list_val/(val_range_max-val_range_min)*(size(vel_color,1)-1)+1);

% Show scatterplot
fig=scatter(vel_abs_img_list_x, vel_abs_img_list_y);

% Scatterplot settings
set(fig,'SizeData', 3); % size of dots
set(fig,'MarkerFacecolor','flat'); % appearance of dots
fig.CData = vel_color(vel_abs_img_list_val,:) % setting colors of individual dots depending on direction
% axis limits
%xlim([xl,xr]);
%ylim([yt,yb]);
% axis labeling
xlabel('Lateral (mm)'); ylabel('Axial (mm)'); % title('Micro-Bubble image');
%set(gca,'Xtick',linspace(0,300,4)); set(gca, 'XTickLabel',linspace(0,3,4));
%set(gca,'Ytick',linspace(0,900,10)); set(gca, 'YTickLabel',linspace(0,9,10));
xlim([1 img_size(2)]);
ylim([1 img_size(1)]);

% save and cofigure current axis handler
set(gca, 'DataAspectRatio',[1 1 1])
set(gca, 'PlotBoxAspectRatio',[1 1 1])
set(gca, 'YDir','reverse')
set(gca, 'Box','on')


% %Colorbar prop
ch = colorbar;
% set(ch,'position',[0.5729 0.5795 0.0306 0.3167]);
 set(ch,'TicksMode','manual');
 set(ch,'Ticks',linspace(0, 1,6));
 set(ch,'TickLabelsMode','manual');
 set(ch,'TickLabels',linspace(val_range_min, val_range_max,6));
 colormap(vel_color);
 ylabel(ch,'Velocity [mm/s]')


%%
% Scatter image from all MB locations
MB_log_copy = MB_log;
MB_index_list = [];

for MB_index = 1:size(MB_log_copy,2)
   if (MB_log_copy(MB_index).age(3) >= MB_age_condition) && (max(MB_log_copy(MB_index).count(:)) <= MB_count_condition)
       MB_index_list = [MB_index_list, sub2ind([img_size(1:2)],MB_log_copy(MB_index).centroid(:,1) ,MB_log_copy(MB_index).centroid(:,2))'];
   end
end

% scatter_matrix = zeros([img_size(1:2)]); 
% for i = 1:size(MB_log_copy,2)
%     scatter_matrix(MB_index_list(i)) = scatter_matrix(MB_index_list(i)) + 1;
% end

[MB_index_y MB_index_x] = ind2sub([img_size(:)]',MB_index_list); 
phantom_index_y = (phantom_positions(:,3)-0.035)/0.03*630;
phantom_index_x = (phantom_positions(:,1)+0.015)/0.03*600;

% Make figure handle
fig_h = figure; clf;
% set position on screen
set(fig_h,'position',[-1850 570 560 420]);

% Show scatterplot with non-zero positions
fig = scatter(MB_index_x, MB_index_y,'r');
% hold on
% scatter(phantom_index_x, phantom_index_y,'b');

% Scatterplot settings
set(fig,'SizeData', 3); % size of dots
set(fig,'MarkerFacecolor','flat'); % appearance of dots

xlabel('Lateral [mm]'); ylabel('Axial [mm]'); % title('Micro-Bubble image');
% set(gca,'Xtick',linspace(0,1189,5)); set(gca, 'XTickLabel',linspace(0,12,5));
% set(gca,'Ytick',linspace(0,2489,6)); set(gca, 'YTickLabel',linspace(0,25,6));
% set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
% set(gca, 'PlotBoxAspectRatio',[1 1 1])
xlim([1 img_size(2)]);
ylim([1 img_size(1)]);

set(gca, 'YDir','reverse'); % reverse y-axis
set(gca, 'Box','on');

