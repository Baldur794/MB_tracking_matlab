%% Calculate movement
% vel_x = mov_x;%----
% vel_y = mov_y;%----
vel_x = vel_x_raw-repmat(mean(vel_x_raw')',1,size(vel_x_raw,2));
vel_y = vel_y_raw-repmat(mean(vel_y_raw')',1,size(vel_y_raw,2));
mov_y = zeros(size(vel_y));
mov_x = zeros(size(vel_x));


% Only for displaying movement compared to ref first image
for idx_wind = 1:size(vel_y,1)
    for idx_frame = start_frame:start_frame+size(mov_y,2)-1;
        mov_y(idx_wind,idx_frame-start_frame+1) = sum(vel_y(idx_wind,1:idx_frame-start_frame+1));
        mov_x(idx_wind,idx_frame-start_frame+1) = sum(vel_x(idx_wind,1:idx_frame-start_frame+1));
    end
end

% figure();plot(abs(fft(vel_y(idx_contrast,:)-mean(vel_y(idx_contrast,:)))));
% figure();plot(abs(fft(vel_x(idx_contrast,:)-mean(vel_x(idx_contrast,:)))));


% Filter out pulse
n = 20;
b_l = fir1(n,0.05,'low');


vel_y_ventilator = filter(b_l,1,vel_y,[],2)*1.22; %1.22 = damping factor from filter
vel_y_ventilator = vel_y_ventilator(:,n/2:end);
vel_x_ventilator = filter(b_l,1,vel_x,[],2)*1.22; %1.22 = damping factor from filter
vel_x_ventilator = vel_x_ventilator(:,n/2:end);

% figure();plot(abs(fft(vel_y_ventilator(idx_contrast,:)-mean(vel_y_ventilator(idx_contrast,:)))));
% figure();plot(abs(fft(vel_x_ventilator(idx_contrast,:)-mean(vel_x_ventilator(idx_contrast,:)))));

% Filter out ventilator
n = 20;
b_h = fir2(n,[0 0.1 0.15 0.3 0.35 1] ,[0 0 1 1 0 0]);
% n = 20;
% b_h = fir1(n,0.05,'high');


vel_y_pulse = filter(b_h,1,vel_y,[],2)*1.2; %1.22 = damping factor from filter
vel_y_pulse = vel_y_pulse(:,n/2:end);
vel_x_pulse = filter(b_h,1,vel_x,[],2)*1.2; %1.22 = damping factor from filter
vel_x_pulse = vel_x_pulse(:,n/2:end);

% figure();plot(abs(fft(vel_y_pulse(idx_contrast,:)-mean(vel_y_pulse(idx_contrast,:)))));
% figure();plot(abs(fft(vel_x_pulse(idx_contrast,:)-mean(vel_x_pulse(idx_contrast,:)))));


%% Interpolate velocity ventilator
interpolate_factor = 10;
f_rep_ventilator = 430;
rep_factor_ventilator = 3%;1000;
vel_yy = spline(1:size(vel_y_ventilator,2),vel_y_ventilator,1:1/interpolate_factor:size(vel_y_ventilator,2));
vel_xx = spline(1:size(vel_x_ventilator,2),vel_x_ventilator,1:1/interpolate_factor:size(vel_x_ventilator,2));

min_variance = 1000;
% for f_rep = 400:450
    % Axial velocity list
    axial_vel_list = zeros(size(vel_y_ventilator,1),f_rep_ventilator,floor(size(vel_yy,2)/f_rep_ventilator));

    for idx_wind = 1:size(vel_yy,1)
        for i = 1:size(axial_vel_list,2)
            for k = 1:size(axial_vel_list,3)
                index = size(axial_vel_list,2) * (k-1) + i;
                axial_vel_list(idx_wind,i,k) = vel_yy(idx_wind,index);
            end
        end
    end

    % Remove mean
    for idx_wind = 1:size(vel_yy,1)
        for k = 1:size(axial_vel_list,3)
            axial_vel_list(idx_wind,:,k) = axial_vel_list(idx_wind,:,k)-mean(axial_vel_list(idx_wind,:,k));
        end
    end

    % Recreated velocity list
    vel_yy_zero_mean = zeros(size(vel_yy));
    for idx_wind = 1:size(vel_yy,1)
        for k = 1:size(axial_vel_list,3)
            vel_yy_zero_mean(idx_wind,size(axial_vel_list,2)*(k-1)+1:size(axial_vel_list,2)*k) = axial_vel_list(idx_wind,:,k);
        end
    end

    % Axial mean
    axial_vel_mean = zeros(size(vel_yy,1),size(axial_vel_list,2));
    axial_vel_mean = mean(axial_vel_list,3);

    % Axial variance
    axial_vel_var = zeros(size(vel_yy,1),1);
    for idx_temp = 1:size(vel_y_ventilator,1)
        axial_vel_var(idx_temp) = sum(var(axial_vel_list(idx_temp,:,:),1,3))/size(axial_vel_list,2);
    end
    
%     if sum(axial_vel_var) < min_variance
%         min_variance = sum(axial_vel_var);
%         f_rep_min = f_rep;
%     end
% end
%
% Axial std
axial_vel_std = zeros(size(vel_yy,1),1);
axial_vel_std = sqrt(axial_vel_var);

axial_vel_mean = repmat(axial_vel_mean,1,rep_factor_ventilator);

% vel_y_mean = spline(linspace(1,f_rep/interpolate_factor,size(axial_vel_mean,2)),axial_vel_mean,1:f_rep/interpolate_factor);
vel_y_mean_ventilator = spline(linspace(1,f_rep_ventilator/interpolate_factor*rep_factor_ventilator,size(axial_vel_mean,2)),axial_vel_mean,1:f_rep_ventilator/interpolate_factor*rep_factor_ventilator);


% Lateral velocity list
lateral_vel_list = zeros(idx_wind,f_rep_ventilator,floor(size(vel_xx,2)/f_rep_ventilator));

for idx_wind = 1:size(vel_yy,1)
    for i = 1:size(lateral_vel_list,2)
        for k = 1:size(lateral_vel_list,3)
            index = size(lateral_vel_list,2) * (k-1) + i;
            lateral_vel_list(idx_wind,i,k) = vel_xx(idx_wind,index);
        end
    end
end

% Remove mean
for idx_wind = 1:size(vel_yy,1)
    for k = 1:size(lateral_vel_list,3)
        lateral_vel_list(idx_wind,:,k) = lateral_vel_list(idx_wind,:,k)-mean(lateral_vel_list(idx_wind,:,k));
    end
end

% Recreated velocity list
vel_xx_zero_mean = zeros(size(vel_xx));
for idx_wind = 1:size(vel_xx,1)
    for k = 1:size(axial_vel_list,3)
        vel_xx_zero_mean(idx_wind,size(axial_vel_list,2)*(k-1)+1:size(axial_vel_list,2)*k) = axial_vel_list(idx_wind,:,k);
    end
end

% lateral mean
lateral_vel_mean = zeros(size(vel_xx,1),size(lateral_vel_list,2));
lateral_vel_mean = mean(lateral_vel_list,3);

% lateral variance
lateral_vel_var = zeros(size(vel_xx,1),1);
for idx_temp = 1:size(vel_xx,1)
    lateral_vel_var(idx_temp) = sum(var(lateral_vel_list(idx_temp,:,:),1,3))/size(lateral_vel_list,2);
end
%
% lateral std
lateral_vel_std = zeros(size(vel_xx,1),1);
lateral_vel_std = sqrt(lateral_vel_var);

lateral_vel_mean = repmat(lateral_vel_mean,1,rep_factor_ventilator);

% vel_x_mean = spline(linspace(1,f_rep/interpolate_factor,size(lateral_vel_mean,2)),lateral_vel_mean,1:f_rep/interpolate_factor);
vel_x_mean_ventilator = spline(linspace(1,f_rep_ventilator/interpolate_factor*rep_factor_ventilator,size(lateral_vel_mean,2)),lateral_vel_mean,1:f_rep_ventilator/interpolate_factor*rep_factor_ventilator);


%% Interpolate velocity pulse
interpolate_factor = 10;
f_rep_pulse = 86;
rep_factor_pulse = rep_factor_ventilator*ceil(f_rep_ventilator/f_rep_pulse);
vel_yy = spline(1:size(vel_y_pulse,2),vel_y_pulse,1:1/interpolate_factor:size(vel_y_pulse,2));
vel_xx = spline(1:size(vel_x_pulse,2),vel_x_pulse,1:1/interpolate_factor:size(vel_x_pulse,2));

% min_variance = 1000;
% for f_rep = 50:150
    % Axial velocity list
    axial_vel_list = zeros(size(vel_y_pulse,1),f_rep_pulse,floor(size(vel_yy,2)/f_rep_pulse));

    for idx_wind = 1:size(vel_yy,1)
        for i = 1:size(axial_vel_list,2)
            for k = 1:size(axial_vel_list,3)
                index = size(axial_vel_list,2) * (k-1) + i;
                axial_vel_list(idx_wind,i,k) = vel_yy(idx_wind,index);
            end
        end
    end

    % Remove mean
    for idx_wind = 1:size(vel_yy,1)
        for k = 1:size(axial_vel_list,3)
            axial_vel_list(idx_wind,:,k) = axial_vel_list(idx_wind,:,k)-mean(axial_vel_list(idx_wind,:,k));
        end
    end

    % Recreated velocity list
    vel_yy_zero_mean = zeros(size(vel_yy));
    for idx_wind = 1:size(vel_yy,1)
        for k = 1:size(axial_vel_list,3)
            vel_yy_zero_mean(idx_wind,size(axial_vel_list,2)*(k-1)+1:size(axial_vel_list,2)*k) = axial_vel_list(idx_wind,:,k);
        end
    end

    % Axial mean
    axial_vel_mean = zeros(size(vel_yy,1),size(axial_vel_list,2));
    axial_vel_mean = mean(axial_vel_list,3);

    % Axial variance
    axial_vel_var = zeros(size(vel_yy,1),1);
    for idx_temp = 1:size(vel_y_pulse,1)
        axial_vel_var(idx_temp) = sum(var(axial_vel_list(idx_temp,:,:),1,3))/size(axial_vel_list,2);
    end
    
%     if sum(axial_vel_var) < min_variance
%         min_variance = sum(axial_vel_var)
%         f_rep_min = f_rep
%     end
% end
%
% Axial std
axial_vel_std = zeros(size(vel_yy,1),1);
axial_vel_std = sqrt(axial_vel_var);

axial_vel_mean = repmat(axial_vel_mean,1,rep_factor_pulse);

% vel_y_mean = spline(linspace(1,f_rep/interpolate_factor,size(axial_vel_mean,2)),axial_vel_mean,1:f_rep/interpolate_factor);
vel_y_mean_pulse = spline(linspace(1,f_rep_pulse/interpolate_factor*rep_factor_pulse,size(axial_vel_mean,2)),axial_vel_mean,1:f_rep_pulse/interpolate_factor*rep_factor_pulse);


% Lateral velocity list
lateral_vel_list = zeros(idx_wind,f_rep_pulse,floor(size(vel_xx,2)/f_rep_pulse));

for idx_wind = 1:size(vel_yy,1)
    for i = 1:size(lateral_vel_list,2)
        for k = 1:size(lateral_vel_list,3)
            index = size(lateral_vel_list,2) * (k-1) + i;
            lateral_vel_list(idx_wind,i,k) = vel_xx(idx_wind,index);
        end
    end
end

% Remove mean
for idx_wind = 1:size(vel_yy,1)
    for k = 1:size(lateral_vel_list,3)
        lateral_vel_list(idx_wind,:,k) = lateral_vel_list(idx_wind,:,k)-mean(lateral_vel_list(idx_wind,:,k));
    end
end

% Recreated velocity list
vel_xx_zero_mean = zeros(size(vel_xx));
for idx_wind = 1:size(vel_xx,1)
    for k = 1:size(axial_vel_list,3)
        vel_xx_zero_mean(idx_wind,size(axial_vel_list,2)*(k-1)+1:size(axial_vel_list,2)*k) = axial_vel_list(idx_wind,:,k);
    end
end


% lateral mean
lateral_vel_mean = zeros(size(vel_xx,1),size(lateral_vel_list,2));
lateral_vel_mean = mean(lateral_vel_list,3);

% lateral variance
lateral_vel_var = zeros(size(vel_xx,1),1);
for idx_temp = 1:size(vel_xx,1)
    lateral_vel_var(idx_temp) = sum(var(lateral_vel_list(idx_temp,:,:),1,3))/size(lateral_vel_list,2);
end
%
% lateral std
lateral_vel_std = zeros(size(vel_xx,1),1);
lateral_vel_std = sqrt(lateral_vel_var);

lateral_vel_mean = repmat(lateral_vel_mean,1,rep_factor_pulse);

% vel_x_mean = spline(linspace(1,f_rep/interpolate_factor,size(lateral_vel_mean,2)),lateral_vel_mean,1:f_rep/interpolate_factor);
vel_x_mean_pulse = spline(linspace(1,f_rep_pulse/interpolate_factor*rep_factor_pulse,size(lateral_vel_mean,2)),lateral_vel_mean,1:f_rep_pulse/interpolate_factor*rep_factor_pulse);




%% Calculate movement (from processed velocity)
mov_y_comp_ventilator = zeros(size(vel_y_mean_ventilator,1),size(vel_y_mean_ventilator,2)+1);
mov_x_comp_ventilator = zeros(size(vel_x_mean_ventilator,1),size(vel_x_mean_ventilator,2)+1);

for idx_wind = 1:size(vel_y_mean_ventilator,1)
    for idx_frame = start_frame:start_frame+size(mov_y_comp_ventilator,2)-2;
        mov_y_comp_ventilator(idx_wind,idx_frame-start_frame+2) = sum(vel_y_mean_ventilator(idx_wind,1:idx_frame-start_frame+1));
        mov_x_comp_ventilator(idx_wind,idx_frame-start_frame+2) = sum(vel_x_mean_ventilator(idx_wind,1:idx_frame-start_frame+1));
    end
end

mov_y_comp_pulse = zeros(size(vel_y_mean_pulse,1),size(vel_y_mean_pulse,2)+1);
mov_x_comp_pulse = zeros(size(vel_x_mean_pulse,1),size(vel_x_mean_pulse,2)+1);

for idx_wind = 1:size(vel_y_mean_pulse,1)
    for idx_frame = start_frame:start_frame+size(mov_y_comp_pulse,2)-2;
        mov_y_comp_pulse(idx_wind,idx_frame-start_frame+2) = sum(vel_y_mean_pulse(idx_wind,1:idx_frame-start_frame+1));
        mov_x_comp_pulse(idx_wind,idx_frame-start_frame+2) = sum(vel_x_mean_pulse(idx_wind,1:idx_frame-start_frame+1));
    end
end

% figure();plot(mov_y_comp(44,:));
% figure();plot(mov_x_comp);

idx_contrast = sub2ind(size(X_comp),4,6);
%% Aligning B-mode and contrast data
offset = 43
% Ventilator
mov_y_comp_contrast_ventilator = [];
mov_x_comp_contrast_ventilator = [];

% Accounting for resolution in B-mode
mov_y_comp_contrast_ventilator = -mov_y_comp_ventilator(idx_contrast,:)*2.55;
mov_x_comp_contrast_ventilator = -mov_x_comp_ventilator(idx_contrast,:)*4.3;

% Aligning data to match contrast data
mov_y_comp_contrast_ventilator = mov_y_comp_contrast_ventilator(1,offset:end);
mov_x_comp_contrast_ventilator = mov_x_comp_contrast_ventilator(1,offset:end);

% resample ventilator
% corr_max = 0;
% for i = 500:1100
%     [P,Q]= rat(i/1000);
%     temp = resample(mov_y_comp_contrast_ventilator,P,Q);
%     if max(crosscorr(axial_mov-mean(axial_mov),temp-mean(temp),100)) > corr_max
%         corr_max = max(crosscorr(axial_mov,temp,100));
%         i_max = i
%     end
% end

[P,Q]= rat(939/1000);
mov_y_comp_contrast_ventilator_resample = resample(mov_y_comp_contrast_ventilator,P,Q);
mov_y_comp_contrast_ventilator_resample = round(mov_y_comp_contrast_ventilator_resample-mean(mov_y_comp_contrast_ventilator_resample));

mov_x_comp_contrast_ventilator_resample = resample(mov_x_comp_contrast_ventilator,P,Q);
mov_x_comp_contrast_ventilator_resample = round(mov_x_comp_contrast_ventilator_resample);

figure(); crosscorr(axial_mov,mov_y_comp_contrast_ventilator_resample,100);
mov_y_comp_contrast_ventilator_resample = -mov_y_comp_contrast_ventilator_resample;
mov_x_comp_contrast_ventilator_resample = -mov_x_comp_contrast_ventilator_resample;

% mov_y_comp_contrast = mov_y_comp_contrast_ventilator_resample;% + mov_y_comp_contrast_pulse(1,1:size(mov_y_comp_contrast_ventilator,2));
% mov_x_comp_contrast = zeros(size(mov_y_comp_contrast_ventilator_resample));%mov_x_comp_contrast_pulse;


%% Pulse
mov_y_comp_contrast_pulse = [];
mov_x_comp_contrast_pulse = [];

mov_y_comp_contrast_pulse = -mov_y_comp_pulse(idx_contrast,:)*2.55;
mov_x_comp_contrast_pulse = -mov_x_comp_pulse(idx_contrast,:)*4.3;

mov_y_comp_contrast_pulse = mov_y_comp_contrast_pulse(1,offset:end);
mov_x_comp_contrast_pulse = mov_x_comp_contrast_pulse(1,offset:end);

% resample pulse ( lateral)
% corr_max = 0;
% for i = 300:700
%     [P,Q]= rat(i/500);
%     temp = resample(mov_x_comp_contrast_pulse,P,Q);
%     if max(crosscorr(lateral_mov,temp,100)) > corr_max
%         corr_max = max(crosscorr(lateral_mov,temp,100));
%         i_max = i;
%     end
% end

% resample pulse (axial)
% corr_max = 0;
% for i = 200:800
%     [P,Q]= rat(i/500);
%     temp = resample(mov_y_comp_contrast_pulse,P,Q);
%     if max(crosscorr(axial_mov,temp,100)) > corr_max
%         corr_max = max(crosscorr(axial_mov,temp,100));
%         i_max = i;
%     end
% end

[P,Q]= rat(797/500);
mov_y_comp_contrast_pulse = resample(mov_y_comp_contrast_pulse,P,Q);
mov_y_comp_contrast_pulse = round(mov_y_comp_contrast_pulse);

mov_x_comp_contrast_pulse = resample(mov_x_comp_contrast_pulse,P,Q);
mov_x_comp_contrast_pulse = round(mov_x_comp_contrast_pulse);

% figure(); crosscorr(lateral_mov,mov_x_comp_contrast_pulse,100);
mov_y_comp_contrast_pulse = -mov_y_comp_contrast_pulse;
mov_x_comp_contrast_pulse = -mov_x_comp_contrast_pulse;

mov_y_comp_contrast = mov_y_comp_contrast_ventilator + mov_y_comp_contrast_pulse(1,1:size(mov_y_comp_contrast_ventilator,2));
mov_x_comp_contrast = mov_x_comp_contrast_ventilator + mov_x_comp_contrast_pulse(1,1:size(mov_x_comp_contrast_ventilator,2));


% figure();plot(mov_y_comp_contrast);
% figure();plot(mov_x_comp_contrast);
% figure();plot(abs(fft(mov_y_comp_contrast-mean(mov_y_comp_contrast))));
% figure();plot(abs(fft(mov_x_comp_contrast-mean(mov_x_comp_contrast))));


%% Displacement axial
figure(); plot((mov_y + flipud(repmat(linspace(0,50,size(mov_y,1))',1,size(mov_y,2))))');% legend('1','2','3','4','5','6','7');
xlabel('frames (50 fps)'); ylabel('Axial displacement (mm)'); 
title('Axial displacement');
ylim([-5,70]);% xlim([0,40]);
%set(gca,'Xtick',linspace(0,100,20)); set(gca,'Xtick',linspace(0,img_size(2),4));
%set(gca,'Ytick',linspace(17,32,5)); set(gca, 'YTickLabel',linspace(-0.2,0.2,5));

%% Velocity axial
figure(); plot((vel_y + flipud(repmat(linspace(0,50,size(vel_y,1))',1,size(vel_y,2))))');% legend('1','2','3','4','5','6','7');
xlabel('frames (50 fps)'); ylabel('Axial Velocity (mm)'); 
title('Axial Velocity');
ylim([-5,70]);% xlim([0,40]);
%set(gca,'Xtick',linspace(0,100,20)); set(gca,'Xtick',linspace(0,img_size(2),4));
%set(gca,'Ytick',linspace(17,32,5)); set(gca, 'YTickLabel',linspace(-0.2,0.2,5));

%% Axial displacement interpolate
figure(); plot((mov_yy + flipud(repmat(linspace(0,40,size(mov_y,1))',1,size(mov_yy,2))))'); legend('1','2','3','4','5','6','7');
%xlabel('frames (50 fps)'); ylabel('Axial displacement (mm)'); 
title('Axial displacement interpolated');
%ylim([-5,45]);% xlim([0,3000]);
%set(gca,'Xtick',linspace(0,100,20)); set(gca,'Xtick',linspace(0,img_size(2),4));
%set(gca,'Ytick',linspace(65.5,95.5,5)); set(gca, 'YTickLabel',linspace(-0.25,0.25,5));

%% Displacement axial mean
figure(); plot((axial_disp_mean)'); legend('1','2','3','4','5','6','7'); 
%xlabel('time (ms)'); ylabel('Axial displacement (mm)');
title('Axial displacement mean');
%ylim([-20,20]);
%set(gca,'Xtick',linspace(0,45,10)); set(gca,'XtickLabel',linspace(0,0.9,10));
%set(gca,'Ytick',linspace(-20,5,5)); set(gca, 'YTickLabel',linspace(-0.16,0.16,5));

%% Displacement multiple axial overlay
for j = 1:size(axial_disp_list,1)
    figure(108);
    
    subplot(size(axial_disp_list,1),1,j)
    
    hold on
    
    for i = 1:size(axial_disp_list,3)
        plot((squeeze(axial_disp_list(j,:,i))));
    end
    hold off
    %legend('1','2','3','4','5','6','7');
    xlabel('Frames'); ylabel('Axial displacement (mm)');
    title('Axial displacement');
    %ylim([-5,5]); %xlim([1 43]);
    %set(gca,'Xtick',linspace(0,37,6)); set(gca,'XtickLabel',linspace(0,20,6));
    %set(gca,'Ytick',linspace(-3,12,5)); set(gca, 'YTickLabel',linspace(0,0.2,5));
end

%% Displacement lateral
figure(); plot((mov_x + flipud(repmat(linspace(0,40,size(mov_x,1))',1,size(mov_x,2))))'); legend('1','2','3','4','5','6','7'); 
xlabel('frames (50 fps)'); ylabel('Lateral displacement (mm)');
title('Lateral displacement');
ylim([-5,55]);% xlim([0,300]);
%set(gca,'Xtick',linspace(0,img_size(2),4)); set(gca,'Xtick',linspace(0,img_size(2),4));
%set(gca,'Ytick',linspace(14,26,5)); set(gca, 'YTickLabel',linspace(-0.25,0.25,5));

%% Velocity lateral
figure(); plot((vel_x + flipud(repmat(linspace(0,40,size(vel_x,1))',1,size(vel_x,2))))'); legend('1','2','3','4','5','6','7'); 
xlabel('frames (50 fps)'); ylabel('Lateral velocity (mm)');
title('Lateral velocity');
ylim([-5,55]);% xlim([0,300]);
%set(gca,'Xtick',linspace(0,img_size(2),4)); set(gca,'Xtick',linspace(0,img_size(2),4));
%set(gca,'Ytick',linspace(14,26,5)); set(gca, 'YTickLabel',linspace(-0.25,0.25,5));

%% Lateral displacement interpolate
figure(); plot((mov_xx + flipud(repmat(linspace(0,40,size(mov_x,1))',1,size(mov_xx,2))))'); legend('1','2','3','4','5','6','7');
%xlabel('frames (50 fps)'); ylabel('Axial displacement (mm)'); 
title('Lateral displacement interpolated');
ylim([-5,45]);% xlim([0,3000]);
%set(gca,'Xtick',linspace(0,100,20)); set(gca,'Xtick',linspace(0,img_size(2),4));
%set(gca,'Ytick',linspace(65.5,95.5,5)); set(gca, 'YTickLabel',linspace(-0.25,0.25,5));


%% Displacement lateral mean
figure(6); plot((lateral_disp_mean)'); legend('1','2','3','4','5','6','7');
%xlabel('time (ms)'); ylabel('Lateral displacement (mm)'); 
title('Lateral displacement mean');
ylim([-2.3,2.3]);
%set(gca,'Xtick',linspace(0,45,10)); set(gca,'XtickLabel',linspace(0,0.9,10));
%set(gca,'Ytick',linspace(-2.3,2.3,5)); set(gca, 'YTickLabel',linspace(-0.1,0.1,5));

%% Displacement multiple lateral overlay
for j = 1:size(axial_disp_list,1)
    figure(107);
    
    subplot(size(axial_disp_list,1),1,j)
    
    hold on
    for i = 1:size(lateral_disp_list,3)
        plot((squeeze(lateral_disp_list(1,:,i))));
    end
    hold off
    %legend('1','2','3','4','5','6','7');
    xlabel('Frames'); ylabel('Lateral displacement (mm)');
    title('Lateral displacement');
    %ylim([-2,2]); %xlim([1 43]);
    %set(gca,'Xtick',linspace(0,37,6)); set(gca,'XtickLabel',linspace(0,20,6));
    %set(gca,'Ytick',linspace(-1.25,1.25,5)); set(gca, 'YTickLabel',linspace(0,0.1,5));
end
%% Frequency spectrum from axial displacement
figure(); plot(linspace(0,fps/2,frames/2),abs(mov_y_frq(:,1:frames/2))'); legend('1','2','3','4','5','6','7');  xlabel('Frequency (Hz)'); ylabel('Magnitude'); % title('Frequency spectrum of axial displacement');
%set(gca,'Xtick',linspace(0,26,25)); set(gca,'Xtick',linspace(0,25,26));
%set(gca,'Ytick',linspace(0,9,10)); set(gca, 'YTickLabel',linspace(0,1,2));

%% Frequency spectrum from axial mean displacement
figure(11); plot(linspace(0,fps/2,size(axial_disp_mean_rep,2)/2),abs(axial_disp_mean_frq(:,1:size(axial_disp_mean_rep,2)/2))'); legend('1','2','3','4','5','6','7');  xlabel('Frequency (Hz)'); ylabel('Magnitude'); % title('Frequency spectrum of mean axial displacement');
%set(gca,'Xtick',linspace(0,26,25)); set(gca,'Xtick',linspace(0,25,26));
%set(gca,'Ytick',linspace(0,9,10)); set(gca, 'YTickLabel',linspace(0,1,2));

%% Frequency spectrum from lateral displacement
figure(); plot(linspace(0,fps/2,frames/2),abs(mov_x_frq(:,1:frames/2))'); legend('1','2','3','4','5','6','7'); xlabel('Frequency (Hz)'); ylabel('Magnitude'); %title('Frequency spectrum of lateral displacement');
%set(gca,'Xtick',linspace(0,26,25)); set(gca,'Xtick',linspace(0,25,26));
%set(gca,'Ytick',linspace(0,9,10)); set(gca, 'YTickLabel',linspace(0,1,2));

%% Frequency spectrum from axial mean displacement
figure(13); plot(linspace(0,fps/2,size(lateral_disp_mean_rep,2)/2),abs(lateral_disp_mean_frq(7,1:size(lateral_disp_mean_rep,2)/2))'); legend('1','2','3','4','5','6','7');  xlabel('Frequency (Hz)'); ylabel('Magnitude'); % title('Frequency spectrum of mean lateral displacement');
%set(gca,'Xtick',linspace(0,26,25)); set(gca,'Xtick',linspace(0,25,26));
%set(gca,'Ytick',linspace(0,9,10)); set(gca, 'YTickLabel',linspace(0,1,2));

%% Video 
img = abs(load_img_B_mode(1));
norm=max(img(:));

for n = 5
    idx_wind = sub2ind(size(X_comp),4,6);
    
    H.img_wind_cord(idx_wind,:) = [400 450 110 140];
    
    h_f=figure;
    outputVideo=VideoWriter(['MB_video_' num2str(n,'%d')]);
    outputVideo.FrameRate=2;
    open(outputVideo);
    
    nframe = 100;
    mov(1:nframe)= struct('cdata',[],'colormap',[]);
    pause(0.2)
    for idx_frame = start_frame:start_frame+nframe
        switch n
            case 1
                % Full img compensated
                img = abs(load_img_B_mode(idx_frame));
                img_comp = imtranslate(img,[round(mov_x(idx_wind,idx_frame-start_frame+1)),round(mov_y(idx_wind,idx_frame-start_frame+1))]);
                limg=20*log10(img_comp/norm);
                hAx = imagesc(limg, [-40 0]); colormap(gray);
                pause(0.2);
            case 2
                % Full img non-compensated
                img = abs(load_img_B_mode(idx_frame));
                limg=20*log10(img/norm);
                hAx = imagesc(limg, [-40 0]); colormap(gray);
                pause(0.2);
            case 3
                % Window compensated avg
                img_new_temp = abs(H.img_new_template(idx_frame, idx_wind, mode));
                img_new_temp_comp = imtranslate(img_new_temp,[(mov_x_comp_contrast_pulse(idx_frame-start_frame+1)),(mov_y_comp_contrast_ventilator(idx_frame-start_frame+1))]);
%                 img_new_temp_comp = imtranslate(img_new_temp,[round(mov_x_comp_contrast(idx_frame-start_frame+1)),round(mov_y_comp_contrast(idx_frame-start_frame+1))]);
                hAx = imagesc(img_new_temp_comp); colormap(gray);
                pause(0.3);
            case 4
                % Window compensated direct
                img_new_temp = abs(H.img_new_template(idx_frame, idx_wind, mode));
                img_new_temp_comp = imtranslate(img_new_temp,[(mov_x(idx_wind,idx_frame-start_frame+1)),(mov_y(idx_wind,idx_frame-start_frame+1))]);
                hAx = imagesc(img_new_temp_comp); colormap(gray);
                pause(0.3);
                
            case 5
                % Window non-compensated
                img_new_temp = abs(H.img_new_template(idx_frame, idx_wind, mode));
                hAx = imagesc(img_new_temp); colormap(gray);
                pause(0.3);
        end
        
        mov=getframe(gcf);
        %mov=getframe(gca);
        writeVideo(outputVideo,mov.cdata);
    end
    close(gcf)
    close(outputVideo);
end

%% Movement to vector list
mov_x_grid = [];
mov_y_grid = [];
for idx_frame = 1:size(mov_y_comp_ventilator,2)
    for idx = 1:size(mov_y_comp_ventilator,1)
        [x,y] = ind2sub(size(X_comp),idx);
        mov_x_grid(idx_frame,x,y) = mov_x_comp_ventilator(idx,idx_frame);
        mov_y_grid(idx_frame,x,y) = mov_y_comp_ventilator(idx,idx_frame);
       
    end
end

vel_x_grid = mov_x_grid;
vel_y_grid = mov_y_grid;
%% Velocity to vector list
vel_x_grid = [];
vel_y_grid = [];
for idx_frame = 1:size(vel_y_mean,2)
    for idx = 1:size(vel_y_mean,1)
        [x,y] = ind2sub(size(X_comp),idx);
        vel_x_grid(idx_frame,x,y) = vel_x_mean(idx,idx_frame);
        vel_y_grid(idx_frame,x,y) = vel_y_mean(idx,idx_frame);
       
    end
end

%% Velocity quiver video
h_f=figure('Position',[-1500 0 1000 1000]);
outputVideo=VideoWriter(['Video_velocity']);
outputVideo.FrameRate=1;
open(outputVideo);
mov(1:200)= struct('cdata',[],'colormap',[]);
pause(0.2)
img = abs(load_img_B_mode(1));
norm = max(img(:));
for i = 1:size(vel_y_grid,1)
    % Full img non-compensated
    img = abs(load_img_B_mode(1));
    limg=20*log10(img/norm);
    hAx = imagesc(limg, [-30 0]); colormap(gray);
    title(['Frame' num2str(i,'%d')]); 
    hold on
    % Quiver plot
    temp_x = squeeze(vel_x_grid(i,:,:));
    temp_y = squeeze(vel_y_grid(i,:,:));
    
    q_plot = quiver(X_comp(:),Y_comp(:),temp_x(:)*5,temp_y(:)*20);
%     q_plot.AutoScaleFactor = 10;
    q_plot.MaxHeadSize = 0.3;
    q_plot.ShowArrowHead = 'on';
    q_plot.AutoScale = 'off';
    q_plot.LineWidth = 2;
    hold off;
    mov=getframe(gcf);
    writeVideo(outputVideo,mov.cdata);
    pause(0.3);
end
close(gcf)
close(outputVideo);
