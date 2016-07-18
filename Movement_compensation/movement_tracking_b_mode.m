%% Displacement RF

% B_mode
img_wind_cord =[700,900,130,150;
                250,450,40,60];
%  img_wind_cord =[200,400,40,60
%                  200,400,125,145
%                  750,950,130,150
%                  1150,1350,30,50
%                  1400,1600,130,150];

% img_wind_cord =[750,850,128,128
%                  750,850,129,129
%                  750,850,130,130
%                  750,850,131,131
%                  750,850,132,132];

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

max_mov_y = 16;
max_mov_x = 3;

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

 % Display Windows

% Dim and corr for boxes
for idx_wind = 1:size(img_wind_cord,1)               
    rectangle_corr(idx_wind,:) = [img_wind_cord(idx_wind,3), img_wind_cord(idx_wind,1), img_wind_cord(idx_wind,4)-img_wind_cord(idx_wind,3), img_wind_cord(idx_wind,2)-img_wind_cord(idx_wind,1)];
end
% Full image
img_disp = load_img_B_mode(1);
figure(1); clf;
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
    text(rectangle_corr(idx_wind,1),rectangle_corr(idx_wind,2)-50, int2str(idx_wind),'Color','r','FontSize',15,'FontWeight','bold');
end
set(gcf,'position',[-1850 570 560 420]);


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
vel_x = [];
vel_y = [];
vel_abs = [];
for idx_wind = 1:size(img_wind_cord,1)
    vel_x(idx_wind,:) = mov_x(idx_wind,2:end)-mov_x(idx_wind,1:end-1);
    vel_y(idx_wind,:) = mov_y(idx_wind,2:end)-mov_y(idx_wind,1:end-1);
    vel_abs(idx_wind,:) = sqrt(vel_x(idx_wind,:).^2+vel_y(idx_wind,:).^2);
end 

%% Calculate fft
mov_y_frq = [];
for idx_wind = 1:size(img_wind_cord,1)
    mov_y_frq(idx_wind,:) = fft(mov_y(idx_wind,:)-mean(mov_y(idx_wind,:)));
end

mov_x_frq = [];
for idx_wind = 1:size(img_wind_cord,1)
    mov_x_frq(idx_wind,:) = fft(mov_x(idx_wind,:)-mean(mov_x(idx_wind,:)));
end

%% Interpolate
interpolate_factor = 11.51;
mov_yy = spline(1:size(mov_y,2),mov_y,1:1/interpolate_factor:size(mov_y,2));
mov_xx = spline(1:size(mov_x,2),mov_x,1:1/interpolate_factor:size(mov_x,2));


% Axial displacement list
temp = [];
for h = 450:510
    f_rep = 493;
    axial_disp_list = zeros(size(mov_y,1),f_rep,floor(size(mov_yy,2)/f_rep));
    
    for idx_wind = 1:size(img_wind_cord,1)
        for i = 1:size(axial_disp_list,2)
            for k = 1:size(axial_disp_list,3)
                index = size(axial_disp_list,2) * (k-1) + i;
                axial_disp_list(idx_wind,i,k) = mov_yy(idx_wind,index);
            end
        end
    end
    
    
    % Axial mean
    axial_disp_mean = zeros(size(img_wind_cord,1),size(axial_disp_list,2));
    axial_disp_mean = sum(axial_disp_list,3)/size(axial_disp_list,3);
    
    % Axial variance
    axial_disp_var = zeros(size(img_wind_cord,1),1);
    for idx_temp = 1:size(img_wind_cord,1)
        axial_disp_var(idx_temp) = sum(var(axial_disp_list(idx_temp,:,:),1,3))/size(axial_disp_list,2);
    end
    %
    % Axial std
    axial_disp_std = zeros(size(img_wind_cord,1),1);
    axial_disp_std = sqrt(axial_disp_var);
    temp(h-449,:) = axial_disp_std ;
end
[t1 t2] = find(temp == min(temp(:)))

mov_y_mean = spline(linspace(1,f_rep/interpolate_factor,size(axial_disp_mean,2)),axial_disp_mean,1:f_rep/interpolate_factor);


% Lateral displacement list
temp = [];
for h = 450:510
    f_rep = 493;
    lateral_disp_list = zeros(idx_wind,f_rep,floor(size(mov_xx,2)/f_rep));
    
    for idx_wind = 1:size(img_wind_cord,1)
        for i = 1:size(lateral_disp_list,2)
            for k = 1:size(lateral_disp_list,3)
                index = size(lateral_disp_list,2) * (k-1) + i;
                lateral_disp_list(idx_wind,i,k) = mov_xx(idx_wind,index);
            end
        end
    end
    
    % lateral mean
    lateral_disp_mean = zeros(size(img_wind_cord,1),size(lateral_disp_list,2));
    lateral_disp_mean = sum(lateral_disp_list,3)/size(lateral_disp_list,3);
    
    % lateral variance
    lateral_disp_var = zeros(size(img_wind_cord,1),1);
    for idx_temp = 1:size(img_wind_cord,1)
        lateral_disp_var(idx_temp) = sum(var(lateral_disp_list(idx_temp,:,:),1,3))/size(lateral_disp_list,2);
    end
    %
    % lateral std
    lateral_disp_std = zeros(size(img_wind_cord,1),1);
    lateral_disp_std = sqrt(lateral_disp_var);
    temp(h-449,:) = lateral_disp_std ;
end
[t3 t4] = find(temp == min(temp(:)))

mov_x_mean = spline(linspace(1,f_rep/interpolate_factor,size(lateral_disp_mean,2)),lateral_disp_mean,1:f_rep/interpolate_factor);
%% 
locs = [];
pks = [];
for idx_wind = 1:size(img_wind_cord,1)
    corr_match = xcorr(mov_y(idx_wind,:),mov_y_mean(idx_wind,:))
    figure(); plot(corr_match);
    findpeaks(corr_match,'NPeaks',floor(size(mov_y,2)/size(mov_y_mean,2)),'SortStr','descend','MinPeakDistance',10);
    % Find peaks
    [pks(idx_wind,:),locs(idx_wind,:)] = findpeaks(corr_match,'NPeaks',floor(size(mov_y,2)/size(mov_y_mean,2)),'SortStr','descend','MinPeakDistance',10);
end
locs = fliplr(locs)

%%
% Make compensation array axial
mov_y_comp = zeros(size(img_wind_cord,1),size(mov_y,2));
for idx_wind = 1:size(img_wind_cord,1)
    for i = 1:size(locs,2)
        mov_y_comp(idx_wind,locs(idx_wind,i)-size(mov_y,2)+1:locs(idx_wind,i)-size(mov_y,2)+size(mov_y_mean,2)) = mov_y_mean(idx_wind,:);
    end
    if size(mov_y_comp,2) > frames
        mov_y_comp(:,frames+1:end) = [];
    end
    figure();
    plot(mov_y_comp(idx_wind,:));
    hold on
    plot(mov_y(idx_wind,:));
    hold off
    xlabel('Frames'); ylabel('Axial displacement (mm)');
end

% figure();
% plot(abs(mov_y_comp(1,:)-mov_y(1,:)));

% Make compensation array lateral
mov_x_comp = zeros(size(img_wind_cord,1),size(mov_x,2));
for idx_wind = 1:size(img_wind_cord,1)
    for i = 1:size(locs,2)
        mov_x_comp(idx_wind,locs(idx_wind,i)-size(mov_x,2)+1:locs(idx_wind,i)-size(mov_x,2)+size(mov_x_mean,2)) = mov_x_mean(idx_wind,:);
    end
    if size(mov_x_comp,2) > frames
        mov_x_comp(:,frames+1:end) = [];
    end
    figure();
    plot(mov_x_comp(idx_wind,:));
    hold on
    plot(mov_x(idx_wind,:));
    hold off
end

% figure();
% plot(abs(mov_x_comp(2,:)-mov_x(2,:)));



%% Displacement axial
figure(); plot((mov_y + flipud(repmat(linspace(0,40,size(mov_y,1))',1,size(mov_y,2))))');% legend('1','2','3','4','5','6','7');
xlabel('frames (50 fps)'); ylabel('Axial displacement (mm)'); 
title('Axial displacement');
ylim([-5,55]);% xlim([0,40]);
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
ylim([-20,20]);
%set(gca,'Xtick',linspace(0,45,10)); set(gca,'XtickLabel',linspace(0,0.9,10));
%set(gca,'Ytick',linspace(-20,5,5)); set(gca, 'YTickLabel',linspace(-0.16,0.16,5));

%% Displacement multiple axial overlay
for j = 1:size(axial_disp_list,1)
    figure(100);
    
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
    figure(101);
    
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

for n = 3:5
    idx_wind = 1;
    h_f=figure;
    outputVideo=VideoWriter(['MB_video_' num2str(n,'%d')]);
    outputVideo.FrameRate=2;
    open(outputVideo);
    
    nframe = 50;
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
                img_new_temp_comp = imtranslate(img_new_temp,[round(mov_x_comp(idx_wind,idx_frame-start_frame+1+49)),round(mov_y_comp(idx_wind,idx_frame-start_frame+1))]);
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
%%
img = abs(load_img_B_mode(1));
norm=max(img(:));
for idx_frame = 1:size(mov_y,2)
    for idx = 1:size(mov_y,1)
        [x,y] = ind2sub(size(X),idx);
        mov_x_grid(idx_frame,x,y) = mov_x(idx,idx_frame);
        mov_y_grid(idx_frame,x,y) = mov_y(idx,idx_frame);
       
    end
end
%%
h_f=figure;
    outputVideo=VideoWriter(['Video']);
    outputVideo.FrameRate=1;
    open(outputVideo);
    mov(1:200)= struct('cdata',[],'colormap',[]);
    pause(0.2)
for i = 1:size(mov_x,2)
    % Full img non-compensated
    img = abs(load_img_B_mode(1));
    limg=20*log10(img/norm);
    hAx = imagesc(limg, [-40 0]); colormap(gray);
    hold on
    % Quiver plot
    q_plot = quiver(X,Y,squeeze(mov_x_grid(i,:,:)),squeeze(mov_y_grid(i,:,:)));
    q_plot.AutoScaleFactor = 0.1;
    q_plot.MaxHeadSize = 0.3;
    q_plot.ShowArrowHead = 'on';
    q_plot.LineWidth = 2;
    hold off;
    mov=getframe(gcf);
        %mov=getframe(gca);
        writeVideo(outputVideo,mov.cdata);
    pause(0.5);
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
LineLength = 0.08;

%some data
[x,y] = meshgrid(0:0.2:2,0:0.2:2);
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