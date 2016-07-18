%% Displacement RF
%  img_wind_cord =[50,100,4,14
%                  140,180,17,27
%                  200,240,10,20
%                  200,240,25,35
%                  90,130,30,40];

%  img_wind_cord =[58,75,6,11
%                  58,76,6,12
%                  58,77,6,13
%                  58,78,6,14
%                  58,79,6,15];

%  img_wind_cord =[118,150,19,24
%                  118,155,19,25
%                  118,160,19,26
%                  118,165,19,27
%                  118,170,19,28];

 img_wind_cord =[70,100,4,12
                 74,105,25,40
                 104,110,26,30
                 133,148,27,29
                 85,95,46,49];

max_mov_y = 4;
max_mov_x = 3;

fps = 47.5;
frames = 250;
start_frame = 7200;
H = block_matching;
H.max_mov_y = max_mov_y;
H.max_mov_x = max_mov_x;
H.cost_function = 'Norm Corr';
H.img_wind_cord = img_wind_cord;
idx_frame = 1;
mode = 'contrast';

% Display Windows

% Dim and corr for boxes
for idx_wind = 1:size(img_wind_cord,1)               
    rectangle_corr(idx_wind,:) = [img_wind_cord(idx_wind,3), img_wind_cord(idx_wind,1), img_wind_cord(idx_wind,4)-img_wind_cord(idx_wind,3), img_wind_cord(idx_wind,2)-img_wind_cord(idx_wind,1)];
end

% Number of frames for avg
n_fore_grnd = frames;
n_bck_grnd = 0;
n_bck_grnd_skip = 1;

img_avg = load_img_contrast(start_frame, n_fore_grnd,n_bck_grnd,n_bck_grnd_skip);

% Remove top of img
img_avg(1:15,:) = 0;

% Apply threshold
threshold = 25;
if threshold ~= 0
    img_avg(img_avg < threshold) = 0;
end

% Structuring element for dilation
se = strel('square',3);
img_avg = imdilate(img_avg,se);


% figure(1); imagesc(img_fix); colormap('gray'); 
figure(2); imagesc(img_avg); colormap('gray'); 
colormap('gray'); xlabel('Lateral (mm)'); ylabel('Axial (mm)'); %title('B-mode image');

% Full image
% img_disp = load_img_contrast(7000,5,1,1);
% figure(2); clf;
% norm = max(abs(img_disp(:)));
% limg=20*log10(abs(img_disp)/norm);
% imagesc(limg,[-30 0]);% xlim([1 size(img,2)]); ylim([1 size(img,1)]);
% colormap('gray'); xlabel('Lateral (mm)'); ylabel('Axial (mm)'); %title('B-mode image');
%set(gca, 'DataAspectRatio',[1 3.46 1]) % set data aspect ratio in zoom box
%set(gca, 'PlotBoxAspectRatio',[1 1 1])
%set(gca,'Xtick',linspace(0,280,5)); set(gca, 'XTickLabel',linspace(0,12,5));
%set(gca,'Ytick',linspace(0,1960,6)); set(gca, 'YTickLabel',linspace(0,25,6));

hold on;
for idx_wind = 1:size(img_wind_cord,1)
    rectangle('position',rectangle_corr(idx_wind,:),'EdgeColor','r','LineWidth', 2);
    text(rectangle_corr(idx_wind,1),rectangle_corr(idx_wind,2)-20, int2str(idx_wind),'Color','r','FontSize',15,'FontWeight','bold');
end
set(gcf,'position',[-1850 120 560 420]);
hold off;


%% Displacement
mov_y = [];
mov_x = [];
threshold = 25;
% Run through all templates
for idx_wind = 1:size(img_wind_cord,1)
    % Load ref img
    idx_frame = start_frame;
    img_ref_wind = H.img_reference_window(idx_frame, idx_wind, img_avg, mode);
    img_ref_wind(img_ref_wind < threshold) = 0;

    % Repeat for all frames
    for idx_frame = start_frame:start_frame+frames-1;        
%         % Load ref img
%         img_ref_wind = H.img_reference_window(idx_frame, idx_wind, img_avg, mode);
%         img_ref_wind(img_ref_wind < threshold) = 0;
        
        % Load new img
        img_new_temp = H.img_new_template(idx_frame, idx_wind, img_avg, mode);
        img_new_temp(img_new_temp < threshold) = 0;
        
        [motion_y motion_x] = H.motion_displacement(img_ref_wind,img_new_temp);
        mov_y(idx_wind,idx_frame-start_frame+1) = motion_y;
        mov_x(idx_wind,idx_frame-start_frame+1) = motion_x;
    end
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

%% Displacement axial
figure(); plot((mov_y + flipud(repmat(linspace(0,40,5)',1,frames)))'); legend('1','2','3','4','5','6','7');
%xlabel('frames (50 fps)'); ylabel('Axial displacement (mm)'); 
title('Axial displacement');
ylim([-5,45]);% xlim([0,40]);
%set(gca,'Xtick',linspace(0,100,20)); set(gca,'Xtick',linspace(0,img_size(2),4));
%set(gca,'Ytick',linspace(65.5,95.5,5)); set(gca, 'YTickLabel',linspace(-0.25,0.25,5));

%% Axial displacement interpolate
figure(); plot((mov_yy + flipud(repmat(linspace(0,40,5)',1,size(mov_yy,2))))'); legend('1','2','3','4','5','6','7');
%xlabel('frames (50 fps)'); ylabel('Axial displacement (mm)'); 
title('Axial displacement interpolated');
ylim([-5,45]);% xlim([0,3000]);
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
figure();
hold on
for i = 1:size(axial_disp_list,3)
    plot((squeeze(axial_disp_list(2,:,i)))); 
end
hold off
%legend('1','2','3','4','5','6','7'); xlabel('time (ms)'); ylabel('Axial displacement (mm)');
title('Axial displacement');
%ylim([-5,5]); %xlim([1 43]);
%set(gca,'Xtick',linspace(0,37,6)); set(gca,'XtickLabel',linspace(0,20,6));
%set(gca,'Ytick',linspace(-14,6,5)); set(gca, 'YTickLabel',linspace(-0.125,0.125,5));

%% Displacement lateral
figure(); plot((mov_x + flipud(repmat(linspace(0,40,5)',1,frames)))'); legend('1','2','3','4','5','6','7'); 
%xlabel('frames (50 fps)'); ylabel('Lateral displacement (mm)');
title('Lateral displacement');
ylim([-5,45]);% xlim([0,300]);
%set(gca,'Xtick',linspace(0,img_size(2),4)); set(gca,'Xtick',linspace(0,img_size(2),4));
%set(gca,'Ytick',linspace(34,46,5)); set(gca, 'YTickLabel',linspace(-0.26,0.26,5));

%% Lateral displacement interpolate
figure(); plot((mov_xx + flipud(repmat(linspace(0,40,5)',1,size(mov_xx,2))))'); legend('1','2','3','4','5','6','7');
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
figure();
hold on
for i = 1:size(lateral_disp_list,3)
    plot((squeeze(lateral_disp_list(2,:,i)))); 
end
hold off
%legend('1','2','3','4','5','6','7'); xlabel('time (ms)'); ylabel('Axial displacement (mm)');
title('Lateral displacement');
ylim([-2,2]); %xlim([1 43]);
%set(gca,'Xtick',linspace(0,37,6)); set(gca,'XtickLabel',linspace(0,20,6));
%set(gca,'Ytick',linspace(-14,6,5)); set(gca, 'YTickLabel',linspace(-0.125,0.125,5));

%% Frequency spectrum from axial displacement
figure(); plot(linspace(0,fps/2,frames/2),abs(mov_y_frq(2,1:frames/2))'); legend('1','2','3','4','5','6','7');  xlabel('Frequency (Hz)'); ylabel('Magnitude'); % title('Frequency spectrum of axial displacement');
%set(gca,'Xtick',linspace(0,26,25)); set(gca,'Xtick',linspace(0,25,26));
%set(gca,'Ytick',linspace(0,9,10)); set(gca, 'YTickLabel',linspace(0,1,2));

%% Frequency spectrum from axial mean displacement
figure(11); plot(linspace(0,fps/2,size(axial_disp_mean_rep,2)/2),abs(axial_disp_mean_frq(:,1:size(axial_disp_mean_rep,2)/2))'); legend('1','2','3','4','5','6','7');  xlabel('Frequency (Hz)'); ylabel('Magnitude'); % title('Frequency spectrum of mean axial displacement');
%set(gca,'Xtick',linspace(0,26,25)); set(gca,'Xtick',linspace(0,25,26));
%set(gca,'Ytick',linspace(0,9,10)); set(gca, 'YTickLabel',linspace(0,1,2));

%% Frequency spectrum from lateral displacement
figure(12); plot(linspace(0,fps/2,frames/2),abs(mov_x_frq(:,1:frames/2))'); legend('1','2','3','4','5','6','7'); xlabel('Frequency (Hz)'); ylabel('Magnitude'); %title('Frequency spectrum of lateral displacement');
%set(gca,'Xtick',linspace(0,26,25)); set(gca,'Xtick',linspace(0,25,26));
%set(gca,'Ytick',linspace(0,9,10)); set(gca, 'YTickLabel',linspace(0,1,2));

%% Frequency spectrum from lateral mean displacement
figure(13); plot(linspace(0,fps/2,size(lateral_disp_mean_rep,2)/2),abs(lateral_disp_mean_frq(7,1:size(lateral_disp_mean_rep,2)/2))'); legend('1','2','3','4','5','6','7');  xlabel('Frequency (Hz)'); ylabel('Magnitude'); % title('Frequency spectrum of mean lateral displacement');
%set(gca,'Xtick',linspace(0,26,25)); set(gca,'Xtick',linspace(0,25,26));
%set(gca,'Ytick',linspace(0,9,10)); set(gca, 'YTickLabel',linspace(0,1,2));



%% Video 
img = abs(load_img_contrast(7000,3,1,1));
norm=max(img(:));

n = 5;
idx_wind = 5;
h_f=figure;
outputVideo=VideoWriter(['MB_video_' mode '_' num2str(n,'%d')]);
outputVideo.FrameRate=2;
open(outputVideo);

nframe = 50;
mov(1:nframe)= struct('cdata',[],'colormap',[]);
pause(0.2)
for idx_frame = start_frame:start_frame+nframe
switch n
    case 1   
        % Full img compensated
        img = abs(load_img_contrast(idx_frame, 3, 1, 1));
        img_comp = imtranslate(img,[-mov_x(idx_wind,idx_frame-start_frame+1),-mov_y(idx_wind,idx_frame-start_frame+1)]);
        hAx = imagesc(img_comp); colormap(gray);
        pause(0.2);   
    case 2
        % Full img non-compensated
        img = abs(load_img_contrast(idx_frame, 3, 1, 1));
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



%%    
threshold = 35;
n_fore_grnd = 1;
n_bck_grnd = 0;
n_bck_grnd_skip = 1;

for i = 7000:7004
    img = load_img_contrast(i,n_fore_grnd,n_bck_grnd,n_bck_grnd_skip);
    if threshold ~= 0
        img(img < threshold) = 0;
    end
    img(img_fix_bw == 0) = 0;
    img(1:15,:) = 0;
    %figure(i); imagesc(img, [30 100]); colormap('gray'); 
    figure(i+1000); plot(img);
    %figure(i+2000); hist(img(img > 0),1000);
end
%% Video

norm=max(img(:));

h_f=figure;
outputVideo=VideoWriter('MB_video_contrast1.avi');
outputVideo.FrameRate=3;
open(outputVideo);

nframe = 50;
nframe_start = 7200;
mov(1:nframe*10)= struct('cdata',[],'colormap',[]);

n_fore_grnd = 10;
n_bck_grnd = 0;
n_bck_grnd_skip = 1;


for i=nframe_start:nframe_start+nframe
    img = load_img_contrast(i, n_fore_grnd, n_bck_grnd, n_bck_grnd_skip);
    limg=20*log10(img/norm);
    %h_b_ax = imagesc(limg, [-30 -20]); colormap(gray);
    h_b_ax = imagesc(img, [20 50]); colormap(gray);
    pause(0.2);
    
    mov=getframe(gcf);
    writeVideo(outputVideo,mov.cdata);
end

close(gcf)
close(outputVideo);

%% Video

n_fore_grnd = 5;
n_bck_grnd = 0;
n_bck_grnd_skip = 10;
threshold = 40;

h_f=figure;
outputVideo=VideoWriter('MB_video_contrast4.avi');
outputVideo.FrameRate=2;
open(outputVideo);

nframe = 100;
nframe_start = 7000;
mov(1:nframe)= struct('cdata',[],'colormap',[]);
pause(0.3);
img_avg = zeros(size(img));
for i=nframe_start:nframe_start+nframe
    img = load_img_contrast(i, n_fore_grnd, n_bck_grnd, n_bck_grnd_skip);
    img(img < threshold) = 0;
%     img(img >= threshold) = 1;
    
    img_avg = img_avg + img;
    h_b_ax = imagesc(img); colormap(gray);
    pause(0.3);

  mov=getframe(gcf);
  %mov=getframe(gca);
  writeVideo(outputVideo,mov.cdata);
end

close(gcf)
close(outputVideo);


%% 
n_bck_grnd = 0;
n_bck_grnd_skip = 10;
threshold = 30;
MB_size = [10,3];
v_MB = 0.5*10^(-3);

nframe = 500;
nframe_start = 4100;

img_avg = zeros(size(img));
for i=nframe_start:nframe_start+nframe
    img = load_img_contrast(i, n_bck_grnd, n_bck_grnd_skip, v_MB);
    img = img(15:end,:);
    img(img < threshold) = 0;
    while any(img(:))
        [img_int img_idx] = max(img(:));
        [y x] = ind2sub(size(img),img_idx);
        img_avg(y,x) = img_avg(y,x) + 1;%img(y,x);
        
        y_start = y-MB_size(1);
        y_end = y+MB_size(1);
        x_start = x-MB_size(2);
        x_end = x-MB_size(2);
        
        if y_start < 1
            y_start = 1;
        end
        if y_end > size(img,1)
            y_end = size(img,1);
        end
        if x_start < 1
            x_start = 1;
        end
        if x_end < size(img,2)
            x_end = size(img,2);
        end
        
        img(y_start:y_end,x_start:x_end) = 0;
        
    end
end

norm = max(img_avg(:));
limg_avg = 20*log10(img_avg/norm);
figure(4); imagesc(limg_avg,[-40 0]); colormap(gray);







