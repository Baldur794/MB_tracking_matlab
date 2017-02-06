for i = 1;%240:100:740
% Use only MBs within min max xy

box_coord_condition_flag = true;
% box_coord_condition = [115, 150, i, i+100]; % [y1, y2, x1, x2]
box_coord_condition = [90 110 190 230]; % [y1, y2, x1, x2] 
MB_box_coord = []; % Filtered list of MBs satisfying conditions

%% Parameters
% Convert from pixel to vel in mm/s

img_resolution = MB_data.imgResolution;
fps = MB_data.fps;

pos_factor_to_mm = img_resolution.lateral_new*1e3;
pos_factor_to_um = img_resolution.lateral_new*1e5;
vel_factor_to_mm = img_resolution.lateral_new*fps*1e3;


% Min distance from regression line
max_distance = 25;
max_velocity = round(100/vel_factor_to_mm); 
p_fit_order_distance = 1;
p_fit_order_velocity = 3;
p_fit_points = 3000;


%% Get pos and vel of MBs
MB_pos_x = [];
MB_pos_y = [];
MB_vel_x = [];
MB_vel_y = [];
for i = 1:size(MB_index_list ,2)
    MB_index = MB_index_list(i);
    MB_pos_x = [MB_pos_x; MB_log_copy(MB_index).old_pos(:,2)];
    MB_pos_y = [MB_pos_y; MB_log_copy(MB_index).old_pos(:,1)];
    MB_vel_x = [MB_vel_x; MB_log_copy(MB_index).vel(:,2)];
    MB_vel_y = [MB_vel_y; MB_log_copy(MB_index).vel(:,1)];
end
% Use only MBs inside defined box
MB_idx = find((MB_pos_y > box_coord_condition(1)) & (MB_pos_y < box_coord_condition(2)) & (MB_pos_x > box_coord_condition(3)) & (MB_pos_x < box_coord_condition(4)));
MB_pos_x = MB_pos_x(MB_idx);
MB_vel_x = MB_vel_x(MB_idx);
MB_pos_y = MB_pos_y(MB_idx);
MB_vel_y = MB_vel_y(MB_idx);

% Calculate absolute velocity
MB_vel_abs = sqrt(MB_vel_x.^2+MB_vel_y.^2);



% Calculate a regrssion line
p_fit = polyfit(MB_pos_x,MB_pos_y,p_fit_order_distance);
x_fit = linspace(min(MB_pos_x),max(MB_pos_x),p_fit_points);
y_fit = polyval(p_fit,x_fit);


%% Flow Profile with all MBs
temp = zeros(size(x_fit,2),1);
dist_list = zeros(size(MB_pos_x,1),1);
non_valid_idx = [];

for i = 1:size(dist_list,1)
    for j = 1:size(temp,1)
        temp(j) = sqrt((MB_pos_x(i)-x_fit(j))^2+(MB_pos_y(i)-y_fit(j))^2);
    end
    % Check if pos is closer than minimum distance
    if min(temp) < max_distance && MB_vel_abs(i) < max_velocity
       
        dist_list(i) = min(temp);
        
        % Calculate line (y=ax+b) for closest point
        min_idx = find(temp == min(temp));
    else
        non_valid_idx = [non_valid_idx i];
    end
end
dist_list(non_valid_idx) = [];
MB_pos_x(non_valid_idx) = [];
MB_pos_y(non_valid_idx) = [];
MB_vel_x(non_valid_idx) = [];
MB_vel_y(non_valid_idx) = [];
MB_vel_abs(non_valid_idx) = [];

% Calculate new regrssion line without far MBs
p_fit = polyfit(MB_pos_x,MB_pos_y,p_fit_order_distance);
x_fit = linspace(min(MB_pos_x),max(MB_pos_x),p_fit_points);
y_fit = polyval(p_fit,x_fit);

%% Plot regrssion line
figure(2);
hold on;
plot(x_fit,y_fit,'--','Linewidth',3);
% xlabel('Lateral [mm]'); ylabel('Axial [mm]'); % title('Micro-Bubble image');
% set(gca,'Xtick',linspace(0,size(lateral_axis,2),6)); set(gca, 'XTickLabel',linspace(round(lateral_axis(1)*1000),round(lateral_axis(end)*1000),6));
% set(gca,'Ytick',linspace(0,size(depth_axis,2),6)); set(gca, 'YTickLabel',linspace(round(depth_axis(1)*1000),round(depth_axis(end)*1000),6));
% set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
% set(gca, 'PlotBoxAspectRatio',[1 1 1])
% set(gca, 'YDir','reverse'); % reverse y-axis
% set(gca, 'Box','on');
%  xlim([0 size(lateral_axis,2)]);
%  ylim([0 size(depth_axis,2)]);

%% Flow Profile with only valid MBs
temp = zeros(size(x_fit,2),1);
dist_list_new = zeros(size(dist_list,1),1);

for i = 1:size(dist_list_new,1)
    for j = 1:size(temp,1)
        temp(j) = sqrt((MB_pos_x(i)-x_fit(j))^2+(MB_pos_y(i)-y_fit(j))^2);
    end
        dist_list_new(i) = min(temp);       
        % Calculate line (y=ax+b) for closest point
        min_idx = find(temp == min(temp));
        if min_idx >= size(temp,1)
            min_idx = size(temp,1)-1;
        end
        a = (y_fit(min_idx+1)-y_fit(min_idx))/(x_fit(min_idx+1)-x_fit(min_idx));
        b = y_fit(min_idx)-a*x_fit(min_idx);
        
        if a*MB_pos_x(i)+b > y_fit(min_idx)
            dist_list_new(i) = -dist_list_new(i);
        end
end

figure, hist(dist_list_new,40);
xlabel('Distance from centrum [um]'); ylabel('Frequency');
set(gca,'Xtick',linspace(-max(abs(dist_list_new)),max(abs(dist_list_new)),5)); set(gca, 'XTickLabel',round(linspace(-max(abs(dist_list_new))*pos_factor_to_um,max(abs(dist_list_new))*pos_factor_to_um,5),0));
% set(gca,'Ytick',linspace(0,max(MB_vel_abs),5)); set(gca, 'YTickLabel',round(linspace(0,max(MB_vel_abs)*vel_factor_to_mm,5),2));
xlim([-max(abs(dist_list_new)) max(abs(dist_list_new))]);

%% Velocity Profile
p_fit_vel = polyfit(dist_list_new,MB_vel_abs,p_fit_order_velocity);
x_fit_vel = linspace(min(dist_list_new),max(dist_list_new),p_fit_points);
y_fit_vel = polyval(p_fit_vel,x_fit_vel);


figure();
plot(dist_list_new,MB_vel_abs,'o','MarkerSize',3,'MarkerEdgeColor','b','MarkerFaceColor','b')
xlabel('Distance from centrum (mm)'); ylabel('Velocity (mm/s)');
set(gca,'Xtick',linspace(-max(abs(dist_list_new)),max(abs(dist_list_new)),5)); set(gca, 'XTickLabel',round(linspace(-max(abs(dist_list_new))*pos_factor_to_mm,max(abs(dist_list_new))*pos_factor_to_mm,5),2));
set(gca,'Ytick',linspace(0,max(MB_vel_abs),5)); set(gca, 'YTickLabel',round(linspace(0,max(MB_vel_abs)*vel_factor_to_mm,5),2));
hold on
plot(x_fit_vel,y_fit_vel,'-r','LineWidth',3);
xlim([-max(abs(dist_list_new)) max(abs(dist_list_new))]);
ylim([10,max_velocity])

%%
% figure;
% hist3([dist_list_new MB_vel_abs],[10,10]);
% set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
% xlabel('Distance from centrum (mm)'); ylabel('Velocity (mm/s)');
% set(gca,'Xtick',linspace(-max(abs(dist_list_new)),max(abs(dist_list_new)),5)); set(gca, 'XTickLabel',round(linspace(-max(abs(dist_list_new))*pos_factor_to_mm,max(abs(dist_list_new))*pos_factor_to_mm,5),2));
% set(gca,'Ytick',linspace(0,max(MB_vel_abs),5)); set(gca, 'YTickLabel',round(linspace(0,max(MB_vel_abs)*vel_factor_to_mm,5),2));


[n,c] = hist3([dist_list_new MB_vel_abs],[10,10]);
n = n';
n_interp = 100;
x_grid = linspace(min(c{1}),max(c{1}),n_interp); %interp1(c{1},c{1},linspace(min(c{1}),max(c{1}),n_interp),'spline');
y_grid = linspace(min(c{2}),max(c{2}),n_interp); %interp1(c{2},c{2},linspace(min(c{2}),max(c{2}),n_interp),'spline');
[X, Y] = meshgrid(c{1},c{2});
[Xq,Yq] = meshgrid(x_grid,y_grid);
z_grid = interp2(X,Y,n,Xq,Yq,'spline');
figure();
surf(x_grid, y_grid, z_grid,'EdgeColor','none');
colormap(hot);
xlabel('Distance from centrum [um]'); ylabel('Velocity [mm/s]'); zlabel('Frequency');
% set(gca,'Xtick',linspace(-max(abs(dist_list_new)),max(abs(dist_list_new)),5)); set(gca, 'XTickLabel',round(linspace(-max(abs(dist_list_new))*pos_factor_to_mm,max(abs(dist_list_new))*pos_factor_to_mm,5),2));
% set(gca,'Ytick',linspace(0,max(MB_vel_abs),5)); set(gca, 'YTickLabel',round(linspace(0,max(MB_vel_abs)*vel_factor_to_mm,5),2));
% xlim([min(x_grid) max(x_grid)]);
% ylim([min(y_grid) max(y_grid)]);
set(gca,'Xtick',linspace(-10,10,5)); set(gca, 'XTickLabel',round(linspace(-100,100,5),2));
set(gca,'Ytick',linspace(0,70/vel_factor_to_mm,8)); set(gca, 'YTickLabel',round(linspace(0,70,8),2));
% xlim([min(x_grid) max(x_grid)]);
% ylim([min(y_grid) max(y_grid)]);
xlim([-10 10]);
ylim([10/vel_factor_to_mm 70/vel_factor_to_mm]);
zlim([min(min(z_grid)) max(max(z_grid))]);

p_fit_vel = polyfit(dist_list_new,MB_vel_abs,p_fit_order_velocity);
x_fit_vel = x_grid;
x_idx = 1:size(x_grid,2);
y_fit_vel = polyval(p_fit_vel,x_fit_vel);
y_idx = zeros(1,size(y_fit_vel,2));
for i = 1:size(y_idx,2)
    [~, idx] = min(abs(y_fit_vel(i)-y_grid));
    y_idx(i) = idx;
end
z_fit_vel = z_grid( sub2ind(size(z_grid), y_idx, x_idx));
% figure();
hold on
plot3(x_fit_vel(x_idx),y_grid(y_idx),z_fit_vel,'--','color','b','LineWidth',2);

end
%% Surf plot

[n,c] = hist3([MB_pos_x MB_pos_y],[10,10]);
n = n';
n_interp = 100;
x_grid = linspace(min(c{1}),max(c{1}),n_interp); %interp1(c{1},c{1},linspace(min(c{1}),max(c{1}),n_interp),'spline');
y_grid = linspace(min(c{2}),max(c{2}),n_interp); %interp1(c{2},c{2},linspace(min(c{2}),max(c{2}),n_interp),'spline');
[X, Y] = meshgrid(c{1},c{2});
[Xq,Yq] = meshgrid(x_grid,y_grid);
z_grid = interp2(X,Y,n,Xq,Yq,'spline');
figure();
surf(x_grid, y_grid, z_grid,'EdgeColor','none');
xlabel('Lateral [mm]'); ylabel('Axial [mm]'); zlabel('Frequency'); % title('Micro-Bubble image');
set(gca,'Xtick',linspace(min(MB_pos_x),max(MB_pos_x),5)); set(gca, 'XTickLabel',round(linspace(min(MB_pos_x)*pos_factor_to_mm,max(MB_pos_x)*pos_factor_to_mm,5),2));
set(gca,'Ytick',linspace(min(MB_pos_y),max(MB_pos_y),5)); set(gca, 'YTickLabel',round(linspace(min(MB_pos_y)*pos_factor_to_mm,max(MB_pos_y)*pos_factor_to_mm,5),2));
set(gca, 'DataAspectRatio',[1 1 5]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])
xlim([min(x_grid) max(x_grid)]);
ylim([min(y_grid) max(y_grid)]);

hold on
h = pcolor(x_grid,y_grid,z_grid);
set(h, 'EdgeColor', 'none');
h.ZData = ones(size(z_grid)) * -max(max(n));
% colormap(hot) % heat map
% colormap(gray);

