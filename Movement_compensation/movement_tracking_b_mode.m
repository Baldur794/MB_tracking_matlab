%% Displacement RF

% B_mode
img_wind_cord = zeros(1,4);
idx_wind = 1;
for i = 1:1:5
    img_wind_cord(idx_wind,:) = [500 900 130+i 130+i];
    idx_wind = idx_wind + 1;
end

% img_wind_cord = zeros(1,4);
% img_wind_cord(1,:) = [500 900 120 140];



% Spinning disk
% 
% img_wind_cord =[180,220,10,50
%                  180,220,350,390
%                  100,140,25,65
%                  10,50,180,220
%                  350,390,180,220];

% img_wind_cord =[370,390,198,198
%                  370,390,199,199
%                  370,390,200,200
%                  370,390,202,202
%                  370,390,203,203];

max_mov_y = 10;
max_mov_x = 3;

fps = 50;
frames = 1000;
start_frame = 18;
H = block_matching;
H.max_mov_y = max_mov_y;
H.max_mov_x = max_mov_x;
H.cost_function = 'xCorr';
H.img_wind_cord = img_wind_cord;
mode = 'b_mode';

 % Display Windows

% Dim and corr for boxes
for idx_wind = 1:size(img_wind_cord,1)               
    rectangle_corr(idx_wind,:) = [img_wind_cord(idx_wind,3), img_wind_cord(idx_wind,1), img_wind_cord(idx_wind,4)-img_wind_cord(idx_wind,3), img_wind_cord(idx_wind,2)-img_wind_cord(idx_wind,1)];
end
% Full image
img_disp = load_img_B_mode(1);
figure(5); clf;
norm = max(abs(img_disp(:)));
limg=20*log10(abs(img_disp)/norm);
imagesc(limg,[-40 0]);% xlim([1 size(img,2)]); ylim([1 size(img,1)]);
colormap('gray'); xlabel('Lateral (mm)'); ylabel('Axial (mm)'); %title('B-mode image');
%set(gca, 'DataAspectRatio',[1 3.46 1]) % set data aspect ratio in zoom box
%set(gca, 'PlotBoxAspectRatio',[1 1 1])
%set(gca,'Xtick',linspace(0,280,5)); set(gca, 'XTickLabel',linspace(0,12,5));
%set(gca,'Ytick',linspace(0,1960,6)); set(gca, 'YTickLabel',linspace(0,25,6));

hold on;
for idx_wind = 1:size(img_wind_cord,1)
    rectangle('position',rectangle_corr(idx_wind,:),'EdgeColor','r','LineWidth', 2);
%     text(rectangle_corr(idx_wind,1),rectangle_corr(idx_wind,2)-50, int2str(idx_wind),'Color','r','FontSize',15,'FontWeight','bold');
end
set(gcf,'position',[-1850 570 560 420]);


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
            % Load ref img
            img_ref_wind = H.img_reference_window(idx_frame, idx_wind, mode);
            
            % Load new img
            img_new_temp = H.img_new_template(idx_frame+1, idx_wind, mode);
            
            [motion_y motion_x] = H.motion_displacement(img_ref_wind,img_new_temp);
            mov_y(idx_wind,idx_frame-start_frame+1) = motion_y;
            mov_x(idx_wind,idx_frame-start_frame+1) = motion_x;
        end        
end

%% Calculate Veloctity (for fixed ref img)
vel_x = [];
vel_y = [];
vel_abs = [];
for idx_wind = 1:size(img_wind_cord,1)
    vel_x(idx_wind,:) = mov_x(idx_wind,2:end)-mov_x(idx_wind,1:end-1);
    vel_y(idx_wind,:) = mov_y(idx_wind,2:end)-mov_y(idx_wind,1:end-1);
    vel_abs(idx_wind,:) = sqrt(vel_x(idx_wind,:).^2+vel_y(idx_wind,:).^2);
end

%% Calculate movement (for variable ref img)
% vel_x = mov_x;%----
% vel_y = mov_y;%----
vel_x = vel_x-repmat(mean(vel_x')',1,size(vel_x,2));
vel_y = vel_y-repmat(mean(vel_y')',1,size(vel_y,2));
mov_y = zeros(size(vel_y));
mov_x = zeros(size(vel_x));

for idx_wind = 1:size(img_wind_cord,1)
    for idx_frame = start_frame:start_frame+size(mov_y,2)-1;
        mov_y(idx_wind,idx_frame-start_frame+1) = sum(vel_y(idx_wind,1:idx_frame-start_frame+1));
        mov_x(idx_wind,idx_frame-start_frame+1) = sum(vel_x(idx_wind,1:idx_frame-start_frame+1));
    end
end

%% Calculate movement (from processed velocity)
mov_y_comp = zeros(size(vel_y_mean,1),size(vel_y_mean,2)+1);
mov_x_comp = zeros(size(vel_x_mean,1),size(vel_x_mean,2)+1);

for idx_wind = 1:size(img_wind_cord,1)
    for idx_frame = start_frame:start_frame+size(mov_y_comp,2)-2;
        mov_y_comp(idx_wind,idx_frame-start_frame+2) = sum(vel_y_mean(idx_wind,1:idx_frame-start_frame+1));
        mov_x_comp(idx_wind,idx_frame-start_frame+2) = sum(vel_x_mean(idx_wind,1:idx_frame-start_frame+1));
    end
end 


%% Calculate movement (from processed velocity, avg)
mov_y_comp_avg = zeros(size(vel_y_mean_avg,1),size(vel_y_mean_avg,2)+1);
mov_x_comp_avg = zeros(size(vel_x_mean_avg,1),size(vel_x_mean_avg,2)+1);

for idx_wind = 1:size(mov_y_comp_avg,1)
    for idx_frame = start_frame:start_frame+size(mov_y_comp_avg,2)-2;
        mov_y_comp_avg(idx_wind,idx_frame-start_frame+2) = sum(vel_y_mean_avg(idx_wind,1:idx_frame-start_frame+1));
        mov_x_comp_avg(idx_wind,idx_frame-start_frame+2) = sum(vel_x_mean_avg(idx_wind,1:idx_frame-start_frame+1));
    end
end 
idx_contrast = sub2ind(size(X),7,5);
mov_y_comp_contrast = [];
mov_x_comp_contrast = [];
mov_y_comp_contrast = mov_y_comp_avg(idx_contrast,:)*51/12.8;
mov_x_comp_contrast = mov_x_comp_avg(idx_contrast,:)*197/43;
mov_y_comp_contrast = mov_y_comp_contrast(1,4:end);
mov_x_comp_contrast = mov_x_comp_contrast(1,4:end);
figure();crosscorr(test,mov_y_comp_contrast,600);

%% Calculate fft
mov_y_frq = [];
for idx_wind = 1:size(img_wind_cord,1)
    mov_y_frq(idx_wind,:) = fft(mov_y(idx_wind,:)-mean(mov_y(idx_wind,:)));
end

mov_x_frq = [];
for idx_wind = 1:size(img_wind_cord,1)
    mov_x_frq(idx_wind,:) = fft(mov_x(idx_wind,:)-mean(mov_x(idx_wind,:)));
end

%% Interpolate velocity
interpolate_factor = 10;
f_rep = 429;
rep_factor = 25;
vel_yy = spline(1:size(vel_y,2),vel_y,1:1/interpolate_factor:size(vel_y,2));
vel_xx = spline(1:size(vel_x,2),vel_x,1:1/interpolate_factor:size(vel_x,2));


% Axial velocity list
axial_vel_list = zeros(size(vel_y,1),f_rep,floor(size(vel_yy,2)/f_rep));

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
axial_vel_mean = zeros(size(img_wind_cord,1),size(axial_vel_list,2));
axial_vel_mean = mean(axial_vel_list,3);

% Axial variance
axial_vel_var = zeros(size(vel_yy,1),1);
for idx_temp = 1:size(vel_y,1)
    axial_vel_var(idx_temp) = sum(var(axial_vel_list(idx_temp,:,:),1,3))/size(axial_vel_list,2);
end
%
% Axial std
axial_vel_std = zeros(size(vel_yy,1),1);
axial_vel_std = sqrt(axial_vel_var);

% Average of nearby lines
axial_vel_mean_avg = zeros(size(axial_vel_mean,1)/rep_x,size(axial_vel_list,2));
for idx_wind = 1:size(axial_vel_mean,1)/rep_x
    axial_vel_mean_avg(idx_wind,:) = mean(axial_vel_mean(rep_x*(idx_wind-1)+1:rep_x*idx_wind,:));
end

axial_vel_mean_avg = repmat(axial_vel_mean_avg,1,rep_factor);

vel_y_mean = spline(linspace(1,f_rep/interpolate_factor,size(axial_vel_mean,2)),axial_vel_mean,1:f_rep/interpolate_factor);
vel_y_mean_avg = spline(linspace(1,f_rep/interpolate_factor*rep_factor,size(axial_vel_mean_avg,2)),axial_vel_mean_avg,1:f_rep/interpolate_factor*rep_factor);


% Lateral velocity list
lateral_vel_list = zeros(idx_wind,f_rep,floor(size(vel_xx,2)/f_rep));

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

% Average of nearby lines
lateral_vel_mean_avg = zeros(size(lateral_vel_mean,1)/rep_x,size(lateral_vel_list,2));
for idx_wind = 1:size(lateral_vel_mean,1)/rep_x
    lateral_vel_mean_avg(idx_wind,:) = mean(lateral_vel_mean(rep_x*(idx_wind-1)+1:rep_x*idx_wind,:));
end

lateral_vel_mean_avg = repmat(lateral_vel_mean_avg,1,rep_factor);

vel_x_mean = spline(linspace(1,f_rep/interpolate_factor,size(lateral_vel_mean,2)),lateral_vel_mean,1:f_rep/interpolate_factor);
vel_x_mean_avg = spline(linspace(1,f_rep/interpolate_factor*rep_factor,size(lateral_vel_mean_avg,2)),lateral_vel_mean_avg,1:f_rep/interpolate_factor*rep_factor);




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

%---
img_wind_cord(1,:) = [720 750 120 140];
H.img_wind_cord(31,:) = [720 750 120 140];
%---

for n = 3:3
    idx_wind = 31;
    h_f=figure;
    outputVideo=VideoWriter(['MB_video_' num2str(n,'%d')]);
    outputVideo.FrameRate=2;
    open(outputVideo);
    
    nframe = 42;
    mov(1:nframe)= struct('cdata',[],'colormap',[]);
    pause(0.2)
    for idx_frame = start_frame+1:start_frame+nframe
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
                img_new_temp_comp = imtranslate(img_new_temp,[round(mov_x_comp(idx_wind,idx_frame-start_frame+1)),round(mov_y_comp(idx_wind,idx_frame-start_frame+1))]);
                hAx = imagesc(img_new_temp_comp); colormap(gray);
                pause(0.3);
            case 4
                % Window compensated direct
                img_new_temp = abs(H.img_new_template(idx_frame, idx_wind, mode));
                img_new_temp_comp = imtranslate(img_new_temp,[round(mov_x(idx_wind,idx_frame-start_frame+1)),round(mov_y(idx_wind,idx_frame-start_frame+1))]);
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
for idx_frame = 1:size(mov_y_mean,2)
    for idx = 1:size(mov_y_mean,1)
        [x,y] = ind2sub(size(X),idx);
        mov_x_grid(idx_frame,x,y) = mov_x_comp(idx,idx_frame);
        mov_y_grid(idx_frame,x,y) = mov_y_comp(idx,idx_frame);
       
    end
end

vel_x_grid = mov_x_grid;
vel_y_grid = mov_y_grid;
%% Velocity to vector list
vel_x_grid = [];
vel_y_grid = [];
for idx_frame = 1:size(vel_y_mean,2)
    for idx = 1:size(vel_y_mean,1)
        [x,y] = ind2sub(size(X),idx);
        vel_x_grid(idx_frame,x,y) = vel_x_mean(idx,idx_frame);
        vel_y_grid(idx_frame,x,y) = vel_y_mean(idx,idx_frame);
       
    end
end

%% Movement to vector list
mov_x_grid = [];
mov_y_grid = [];
for idx_frame = 1:size(mov_y_comp_avg,2)
    for idx = 1:size(mov_y_comp_avg,1)
        [x,y] = ind2sub(size(X),idx);
        mov_x_grid(idx_frame,x,y) = mov_x_comp_avg(idx,idx_frame);
        mov_y_grid(idx_frame,x,y) = mov_y_comp_avg(idx,idx_frame);
       
    end
end

vel_x_grid = mov_x_grid;
vel_y_grid = mov_y_grid;
%% Velocity to vector list
vel_x_grid = [];
vel_y_grid = [];
for idx_frame = 1:size(vel_y_mean_avg,2)
    for idx = 1:size(vel_y_mean_avg,1)
        [x,y] = ind2sub(size(X),idx);
        vel_x_grid(idx_frame,x,y) = vel_x_mean_avg(idx,idx_frame);
        vel_y_grid(idx_frame,x,y) = vel_y_mean_avg(idx,idx_frame);
       
    end
end

%% Velocity quiver video
h_f=figure('Position',[-1500 0 1000 1000]);
outputVideo=VideoWriter(['Video_velocity']);
outputVideo.FrameRate=1;
open(outputVideo);
mov(1:200)= struct('cdata',[],'colormap',[]);
pause(0.2)
for i = 1:size(vel_y_grid,1)
    % Full img non-compensated
    img = abs(load_img_B_mode(1));
    limg=20*log10(img/norm);
    hAx = imagesc(limg, [-40 0]); colormap(gray);
    title(['Frame' num2str(i,'%d')]); 
    hold on
    % Quiver plot
    q_plot = quiver(X,Y,squeeze(vel_x_grid(i,:,:))*5,squeeze(vel_y_grid(i,:,:))*20);
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


%%

headWidth = 5;
headLength = 5;
LineLength = 0.05;

figure(102); clf
xlim([0 3])
ylim([0 2])
%set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
%set(gca, 'PlotBoxAspectRatio',[1 1 1])
%axis off;
title('Quiver - annotations ','FontSize',16);
%hold on;
%for i = 1:size(mov_x,1)
    for ii = 1:size(X,1)
        for ij = 1:size(X,2)                  
            arrow_annotation = annotation('arrow','headStyle','plain','HeadLength',headLength,'HeadWidth',headWidth);
            set(arrow_annotation,'parent',gca);
            set(arrow_annotation,'position',[X(ii,ij)/1000 Y(ii,ij)/1000 LineLength*mov_x(i,ii,ij) LineLength*mov_y(i,ii,ij)]);
            
        end
    end
%     pause(0.5)
%     clf
%     xlim([1 280])
%     ylim([1 1960])
%     set(gca, 'DataAspectRatio',[1 4 1]) % set data aspect ratio in zoom box
%     set(gca, 'PlotBoxAspectRatio',[1 1 1])
%end

%%

headWidth = 8;
headLength = 8;
LineLength = 0.8;

%some data
[x,y] = meshgrid(0:20:200,0:20:200);
u = cos(x).*y;
v = sin(x).*y;

%quiver plots
figure('Position',[10 10 1000 600],'Color','w');
hax_1 = subplot(1,2,1);
hq = quiver(x,y,u,v);           %get the handle of quiver
title('Regular Quiver plot','FontSize',16);

%get the data from regular quiver
U = hq.UData;
V = hq.VData;
X = hq.XData;
Y = hq.YData;

%right version (with annotation)
hax_2 = subplot(1,2,2);
%hold on;
for ii = 1:length(X)
    for ij = 1:length(X)

        headWidth = 5;
        ah = annotation('arrow',...
            'headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
        set(ah,'parent',gca);
        set(ah,'position',[X(ii,ij) Y(ii,ij) LineLength*U(ii,ij) LineLength*V(ii,ij)]);

    end
end
%axis off;
title('Quiver - annotations ','FontSize',16);

linkaxes([hax_1 hax_2],'xy');

%%
fig_h = figure();

img = abs(load_img_B_mode(1));
img_h = imagesc(img);
ax_h = gca;
ax_h.YDir = 'normal';
%ax_h.PlotBoxAspectRatio = [1 1 1]
ax_h.DataAspectRatio = [250 250 1]

arrow_annotation = annotation('arrow')
set(arrow_annotation,'parent',gca);
arrow_annotation.Color = 'r';
arrow_annotation.Position = [100 100 10 10];

%ax_h.YDir = 'reverse';


% %% Interpolate movement
% interpolate_factor = 11.51;
% mov_yy = spline(1:size(mov_y,2),mov_y,1:1/interpolate_factor:size(mov_y,2));
% mov_xx = spline(1:size(mov_x,2),mov_x,1:1/interpolate_factor:size(mov_x,2));
% 
% 
% % Axial displacement list
% temp = [];
% for h = 450:510
%     f_rep = 493;
%     axial_disp_list = zeros(size(mov_y,1),f_rep,floor(size(mov_yy,2)/f_rep));
%     
%     for idx_wind = 1:size(mov_yy,1)
%         for i = 1:size(axial_disp_list,2)
%             for k = 1:size(axial_disp_list,3)
%                 index = size(axial_disp_list,2) * (k-1) + i;
%                 axial_disp_list(idx_wind,i,k) = mov_yy(idx_wind,index);
%             end
%         end
%     end
%     
%     % Remove mean
%     for idx_wind = 1:size(mov_yy,1)
%         mean_temp = mean(mean(axial_disp_list(idx_wind,:,:)));
%         for k = 1:size(axial_disp_list,3)
%             axial_disp_list(idx_wind,:,k) = axial_disp_list(idx_wind,:,k)-mean(axial_disp_list(idx_wind,:,k));
%         end
%     end
%     
%     % Axial mean
%     axial_disp_mean = zeros(size(img_wind_cord,1),size(axial_disp_list,2));
%     axial_disp_mean = mean(axial_disp_list,3);
%     
%     % Axial variance
%     axial_disp_var = zeros(size(mov_yy,1),1);
%     for idx_temp = 1:size(mov_y,1)
%         axial_disp_var(idx_temp) = sum(var(axial_disp_list(idx_temp,:,:),1,3))/size(axial_disp_list,2);
%     end
%     %
%     % Axial std
%     axial_disp_std = zeros(size(mov_yy,1),1);
%     axial_disp_std = sqrt(axial_disp_var);
%     temp(h-449,:) = axial_disp_std ;
% end
% [t1 t2] = find(temp == min(temp(:)));
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
%     for idx_wind = 1:size(mov_yy,1)
%         for i = 1:size(lateral_disp_list,2)
%             for k = 1:size(lateral_disp_list,3)
%                 index = size(lateral_disp_list,2) * (k-1) + i;
%                 lateral_disp_list(idx_wind,i,k) = mov_xx(idx_wind,index);
%             end
%         end
%     end
%     
%     % Remove mean
%     for idx_wind = 1:size(mov_yy,1)
%         mean_temp = mean(mean(lateral_disp_list(idx_wind,:,:)));
%         for k = 1:size(lateral_disp_list,3)
%             lateral_disp_list(idx_wind,:,k) = lateral_disp_list(idx_wind,:,k)-mean(lateral_disp_list(idx_wind,:,k));
%         end
%     end
%     
%     % lateral mean
%     lateral_disp_mean = zeros(size(mov_yy,1),size(lateral_disp_list,2));
%     lateral_disp_mean = mean(lateral_disp_list,3);
%     
%     % lateral variance
%     lateral_disp_var = zeros(size(mov_yy,1),1);
%     for idx_temp = 1:size(mov_yy,1)
%         lateral_disp_var(idx_temp) = sum(var(lateral_disp_list(idx_temp,:,:),1,3))/size(lateral_disp_list,2);
%     end
%     %
%     % lateral std
%     lateral_disp_std = zeros(size(mov_yy,1),1);
%     lateral_disp_std = sqrt(lateral_disp_var);
%     temp(h-449,:) = lateral_disp_std ;
% end
% [t3 t4] = find(temp == min(temp(:)));
% 
% mov_x_mean = spline(linspace(1,f_rep/interpolate_factor,size(lateral_disp_mean,2)),lateral_disp_mean,1:f_rep/interpolate_factor);

% %% Repeatead mean
% mov_y_comp = zeros(size(mov_yy,1),size(axial_disp_list,2)*size(axial_disp_list,3));
% for idx_wind = 1:size(mov_yy,1)
%     for i = 1:size(axial_disp_list,3)
%         mov_y_comp(idx_wind,size(axial_disp_list,2)*(i-1)+1:size(axial_disp_list,2)*i) = axial_disp_mean(idx_wind,:); 
%     end
% end
% mov_y_comp = spline(linspace(1,f_rep/interpolate_factor*size(axial_disp_list,3),size(mov_y_comp,2)),mov_y_comp,1:f_rep/interpolate_factor*size(axial_disp_list,3));
% 
% mov_x_comp = zeros(size(mov_xx,1),size(lateral_disp_list,2)*size(lateral_disp_list,3));
% for idx_wind = 1:size(mov_xx,1)
%     for i = 1:size(lateral_disp_list,3)
%         mov_x_comp(idx_wind,size(lateral_disp_list,2)*(i-1)+1:size(lateral_disp_list,2)*i) = lateral_disp_mean(idx_wind,:); 
%     end
% end
% mov_x_comp = spline(linspace(1,f_rep/interpolate_factor*size(lateral_disp_list,3),size(mov_x_comp,2)),mov_x_comp,1:f_rep/interpolate_factor*size(lateral_disp_list,3));

% %% Match by correlation 
% locs = [];
% pks = [];
% for idx_wind = 1:size(mov_yy,1)
%     corr_match = xcorr(mov_y(idx_wind,:),mov_y_mean(idx_wind,:))
%     figure(); plot(corr_match);
%     findpeaks(corr_match,'NPeaks',floor(size(mov_y,2)/size(mov_y_mean,2)),'SortStr','descend','MinPeakDistance',10);
%     % Find peaks
%     [pks(idx_wind,:),locs(idx_wind,:)] = findpeaks(corr_match,'NPeaks',floor(size(mov_y,2)/size(mov_y_mean,2)),'SortStr','descend','MinPeakDistance',10);
% end
% locs = fliplr(locs)

%%
% % Make compensation array axial
% mov_y_comp = zeros(size(mov_yy,1),size(mov_y,2));
% for idx_wind = 1:size(mov_yy,1)
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
% mov_x_comp = zeros(size(mov_yy,1),size(mov_x,2));
% for idx_wind = 1:size(mov_yy,1)
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
% 
% % figure();
% % plot(abs(mov_x_comp(2,:)-mov_x(2,:)));
% 

%% Displacement multiple axial overlay mean
  figure();
for j = 1:size(mov_y_comp,1)
    
    subplot(size(mov_y_comp,1),1,j)
    
    hold on
    
    for i = 1:size(mov_y_comp,3)
        plot((squeeze(mov_y_comp(j,:,i))));
    end
    hold off
    ylim([-2,8]); %xlim([1 43]);
    set(gca,'Xtick',linspace(1,40,5)); set(gca,'XtickLabel',linspace(0,0.8,5));
    set(gca,'Ytick',linspace(-2,8,5)); set(gca, 'YTickLabel',linspace(-64,64,5));
end
xlabel('Time (s)');
subplot(size(mov_y_comp,1),1,3)
ylabel('Axial displacement (\mum)');
subplot(size(mov_y_comp,1),1,1)
title('Axial displacement mean');

%% Displacement multiple lateral overlay mean
  figure();
for j = 1:size(mov_x_comp,1)
    
    subplot(size(mov_x_comp,1),1,j)
    
    hold on
    
    for i = 1:size(mov_x_comp,3)
        plot((squeeze(mov_x_comp(j,:,i))));
    end
    hold off
    ylim([-2,1]); %xlim([1 43]);
    set(gca,'Xtick',linspace(1,40,5)); set(gca,'XtickLabel',linspace(0,0.8,5));
    set(gca,'Ytick',linspace(-2,1,5)); set(gca, 'YTickLabel',linspace(-64,64,5));
end
xlabel('Time (s)');
subplot(size(mov_x_copy,1),1,3)
ylabel('Lateral displacement (\mum)');
subplot(size(mov_x_copy,1),1,1)
title('Lateral displacement mean');


%%
%% Displacement multiple axial overlay mean
  figure();
for j = 1:5
    
    subplot(5,1,j)
    
    hold on
    
    for i = 1:size(mov_y_comp,3)
        plot((squeeze(mov_y_comp(j,:,i))));
    end
    hold off
    %ylim([-2,8]); %xlim([1 43]);
    %set(gca,'Xtick',linspace(1,40,5)); set(gca,'XtickLabel',linspace(0,0.8,5));
    %set(gca,'Ytick',linspace(-2,8,5)); set(gca, 'YTickLabel',linspace(-64,64,5));
end

%% Displacement multiple axial overlay mean
  figure();
for j = 1:5
    
    subplot(5,1,j)
    
    hold on
    
    for i = 1:size(mov_y_comp_avg,3)
        plot((squeeze(mov_y_comp_avg(j,:,i))));
    end
    hold off
    %ylim([-2,8]); %xlim([1 43]);
    %set(gca,'Xtick',linspace(1,40,5)); set(gca,'XtickLabel',linspace(0,0.8,5));
    %set(gca,'Ytick',linspace(-2,8,5)); set(gca, 'YTickLabel',linspace(-64,64,5));
end
%xlabel('Time (s)');
%subplot(size(mov_y_copy,1),1,3)
%ylabel('Axial displacement (\mum)');
%subplot(size(mov_y_copy,1),1,1)
%title('Axial displacement mean');
