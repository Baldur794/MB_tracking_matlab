vel_abs_img_list_x_temp = vel_abs_img_list_x;
vel_abs_img_list_y_temp = vel_abs_img_list_y;
vel_abs_img_list_val_temp = vel_abs_img(vel_abs_img_list);

p_fit = polyfit(vel_abs_img_list_x,vel_abs_img_list_y,2);
x_fit = linspace(64,260,3000);
y_fit = polyval(p_fit,x_fit);
%%
vel_abs_img_list_x_temp = vel_abs_img_list_x;
vel_abs_img_list_y_temp = vel_abs_img_list_y;
vel_abs_img_list_val_temp = vel_abs_img(vel_abs_img_list);


dist_list_temp = zeros(size(x_fit,2),1);
dist_list = zeros(size(vel_abs_img_list_x_temp,1),1);
% for i = 1:size(dist_list,1)
%     for j = 1:size(dist_list_temp,1)
%         dist_list_temp(j) = sqrt((vel_abs_img_list_y(i)-y_fit(j))^2+(vel_abs_img_list_x(i)-x_fit(j))^2);
%     end
%     dist_list(i) = min(dist_list_temp);
% end
% figure, hist(dist_list,20)

for i = 1:size(dist_list,1)
    for j = 1:size(dist_list_temp,1)
        dist_list_temp(j) = sqrt((vel_abs_img_list_y(i)-y_fit(j))^2+(vel_abs_img_list_x(i)-x_fit(j))^2);
    end
    if min(dist_list_temp) < 20
        dist_list(i) = min(dist_list_temp);
    else
        dist_list(i) = -1;
    end
end
vel_abs_img_list_x_temp(find(dist_list == -1)) = [];
vel_abs_img_list_y_temp(find(dist_list == -1)) = [];
vel_abs_img_list_val_temp(find(dist_list == -1)) = [];
dist_list(find(dist_list == -1)) = [];
p_fit = polyfit(vel_abs_img_list_x_temp,vel_abs_img_list_y_temp,3);
x_fit = linspace(64,260,1000);
y_fit = polyval(p_fit,x_fit);

%figure, hist(dist_list,20)


%%
dist_list_temp = zeros(size(x_fit,2),1);
dist_list = zeros(size(vel_abs_img_list_x_temp,1),1);
for i = 1:size(dist_list,1)
    for j = 1:size(dist_list_temp,1)
        dist_list_temp(j) = sqrt((vel_abs_img_list_y_temp(i)-y_fit(j))^2+(vel_abs_img_list_x_temp(i)-x_fit(j))^2);
        if (vel_abs_img_list_x_temp(i)-x_fit(j)) < 0
            dist_list_temp(j) = -dist_list_temp(j);
        end
    end
    [dist_list(i) dist_list_index] = min(abs(dist_list_temp));
    dist_list(i) = dist_list(i)*sign(dist_list_temp(dist_list_index));
end
figure, hist(dist_list,40)
xlabel('Distance from centrum (mm)'); ylabel('Number of micro-bubbles');
set(gca,'Xtick',linspace(-20,20,5)); set(gca, 'XTickLabel',linspace(-0.2,0.2,5));
xlim([-20 20]);


%%
p_fit_vel = polyfit(dist_list,vel_abs_img_list_val_temp,2);
x_fit_vel = linspace(-9,9,300);
y_fit_vel = polyval(p_fit_vel,x_fit_vel);

figure();
plot(dist_list,vel_abs_img_list_val_temp,'o','MarkerSize',3,'MarkerEdgeColor','b','MarkerFaceColor','b')
xlabel('Distance from centrum (mm)'); ylabel('Velocity (mm/s)');
set(gca,'Xtick',linspace(-10,10,5)); set(gca, 'XTickLabel',linspace(-0.1,0.1,5));
hold on
plot(x_fit_vel,y_fit_vel,'-r','LineWidth',2);
xlim([-10 10]);
ylim([0,40])
%%
hold on
plot(x_fit,y_fit,'-b','LineWidth',2);

%%
