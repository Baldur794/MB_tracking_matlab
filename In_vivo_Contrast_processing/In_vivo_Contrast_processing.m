%% 
MB_window_coord_search = zeros(2,2); % Coordinates for search window
MB_window_coord_search_clear = zeros(2,2); % Coordinates for search clear window

MB_window_size_search_localization = [30 30]; % [y,x] Search window for localization of PSF's
MB_window_size_search_new = [7 7]; % [y,x] Search window for new MB's
MB_window_size_search_existing = [5 5]; % [y,x] Search window for ''old'' MB's
MB_window_size_search_clear = [20 20]; % [y,x] Search window for ''old'' MB's


MB_age_condition = 1;
MB_count_condition = 40; % Number of areas within search windows
MB_window_size_density_avg = [5 3]; % Window for density condition
MB_window_size_density_stuck = [10 10]; % Window for density condition
MB_window_size_density_single = [5 5]; % Window for density condition
MB_dens_condition_avg = 0; % Density condition average
MB_dens_condition_stuck = 40;
MB_dens_condition_single = 0; % Density condition single
MB_dens_condition_single_avg = 0;
MB_vel_mean_condition_x = 0;
MB_vel_mean_condition_y = 0;
MB_vel_variance_condition = 20; % Dont touch...
MB_window_threshold = 2; % Window Threshold (actual MB_window_threshold = Max_intensity*1/threshold)
weighing_factor = 1; % Distance weighing factor
weighing_filter_radius = 3; % Radius around centroid to be considered

MB_window_out_of_bounce = 0; % Checks if search windows is outside image

fps = 47;
n_bck_grnd = 10; 
n_bck_grnd_skip =30;
v_MB = 0.5*10^(-3);
nframe = 100;
idx_frame_start = 4000;
img_size = [2489,1183];
interpolation_type = 'spline';


%%
% initializing MB struct
MB_index = 0; % Index for each MB
MB = []; % Current MB
MB.state = 0;
MB.old_pos = zeros(1,2);
MB.new_pos = zeros(1,2);
MB.vel = zeros(1,2); 
MB.max_int = 0;
MB.centroid = zeros(1,2);
MB.id = 0;
MB.count = 0;
MB_log = MB;


img = load_img_contrast(1000,n_bck_grnd,n_bck_grnd_skip,0.5*10^-3);
[X,Y] = meshgrid(1:size(img,2),1:size(img,1));
[Xq,Yq] = meshgrid(1:1/19.7:size(img,2),1:1/5.1:size(img,1));

tic
for idx_frame=idx_frame_start:idx_frame_start+nframe
    idx_frame
    % Load img
    img = load_img_contrast(idx_frame,n_bck_grnd,n_bck_grnd_skip,v_MB);
    img = interp2(X,Y,img,Xq,Yq,interpolation_type);
%     figure(); imagesc(img); colormap('gray');%---
    
    % Calculate threshold
    SEM = std(img(:));               
    ts = tinv(0.9999999999,length(img(:))-1);     
    CI = mean(img(:)) + ts*SEM;
    global_threshold = CI;
    
    % Threshold
    img_global_threshold = img;
    img_global_threshold(img_global_threshold < global_threshold) = 0;
    %figure(); imagesc(img_global_threshold); colormap('gray');%---
    
    % blob labeling from threshold img
    [img_blob_label_global,blob_count_global] = bwlabel(img_global_threshold,4);
    
    % Sort MB's after intensity
    MB_sort_index = logical([MB.state]);
    [~,MB_sort_index] = sort([MB.max_int].*MB_sort_index,'descend'); 
    % Connect MB's
    for i = 1:sum([MB.state])
        % Update MB_index
        MB_index = MB_sort_index(i);

        % Update window coordinates
         % Update window coordinates
        if MB(MB_index).age(3) == 1
            MB_window_coord_search(1,1) = MB(MB_index).new_pos(1)-MB_window_size_search_new(1);
            MB_window_coord_search(2,1) = MB(MB_index).new_pos(1)+MB_window_size_search_new(1);
            MB_window_coord_search(1,2) = MB(MB_index).new_pos(2)-MB_window_size_search_new(2);
            MB_window_coord_search(2,2) = MB(MB_index).new_pos(2)+MB_window_size_search_new(2);
        else
            MB_window_coord_search(1,1) = MB(MB_index).new_pos(1)-MB_window_size_search_existing(1)+MB(MB_index).vel(1);
            MB_window_coord_search(2,1) = MB(MB_index).new_pos(1)+MB_window_size_search_existing(1)+MB(MB_index).vel(1);
            MB_window_coord_search(1,2) = MB(MB_index).new_pos(2)-MB_window_size_search_existing(2)+MB(MB_index).vel(2);
            MB_window_coord_search(2,2) = MB(MB_index).new_pos(2)+MB_window_size_search_existing(2)+MB(MB_index).vel(2);
        end
        
        % Check for out of bounce
        % y_start
        MB_window_out_of_bounce = 0;
        if MB_window_coord_search(1,1) <= 0
            MB_window_coord_search(1,1) = 1;
            MB_window_out_of_bounce = 1;
        end
        % y_end
        if MB_window_coord_search(2,1) > size(img,1)
            MB_window_coord_search(2,1) = size(img,1);
            MB_window_out_of_bounce = 1;
        end
        % x_start
        if MB_window_coord_search(1,2) <= 0
            MB_window_coord_search(1,2) = 1;
            MB_window_out_of_bounce = 1;
        end
        % x_end
        if MB_window_coord_search(2,2) > size(img,2)
            MB_window_coord_search(2,2) = size(img,2);
            MB_window_out_of_bounce = 1;
        end
        
        % Update clear window coordinates
        MB_window_coord_search_clear(1,1) = MB(MB_index).new_pos(1)-MB_window_size_search_clear(1);
        MB_window_coord_search_clear(2,1) = MB(MB_index).new_pos(1)+MB_window_size_search_clear(1);
        MB_window_coord_search_clear(1,2) = MB(MB_index).new_pos(2)-MB_window_size_search_clear(2);
        MB_window_coord_search_clear(2,2) = MB(MB_index).new_pos(2)+MB_window_size_search_clear(2);
        
        % Check for out of bounce
        % y_start      
        if MB_window_coord_search_clear(1,1) <= 0
            MB_window_coord_search_clear(1,1) = 1;
        end
        % y_end
        if MB_window_coord_search_clear(2,1) > size(img,1)
            MB_window_coord_search_clear(2,1) = size(img,1);
        end
        % x_start
        if MB_window_coord_search_clear(1,2) <= 0
            MB_window_coord_search_clear(1,2) = 1;
        end
        % x_end
        if MB_window_coord_search_clear(2,2) > size(img,2)
            MB_window_coord_search_clear(2,2) = size(img,2);
        end
        
        % Create temp window
        img_temp_window = img_global_threshold(MB_window_coord_search(1,1):MB_window_coord_search(2,1),MB_window_coord_search(1,2):MB_window_coord_search(2,2));
        %figure(); imagesc(img_temp_window); colormap('gray');%--- 
        
         % Check if any blobs are within window
        if any(img_temp_window(:))
            
            % Create temp window
            img_temp_window = img(MB_window_coord_search(1,1):MB_window_coord_search(2,1),MB_window_coord_search(1,2):MB_window_coord_search(2,2));
            %figure(); imagesc(img_temp_window); colormap('gray');%---
            
            % Find max intensity within window
            [max_int, max_index] = max(img_temp_window(:));
            [max_y_window, max_x_window] = ind2sub(size(img_temp_window), max_index);
            max_y = max_y_window + MB_window_coord_search(1,1)-1;
            max_x = max_x_window + MB_window_coord_search(1,2)-1;          

            % Thresholding
            img_temp_window(find(img_temp_window < max_int/MB_window_threshold)) = 0;
            %figure(); imagesc(img_temp_window); colormap('gray');%---
            
            % Count number of blobs in window from global img
            img_blob_label_global_window = img_blob_label_global(MB_window_coord_search_clear(1,1):MB_window_coord_search_clear(2,1),MB_window_coord_search_clear(1,2):MB_window_coord_search_clear(2,2));
            blobs_global_window = unique(img_blob_label_global_window);
            blobs_global_window(1) = [];
            blob_count_global_window = length(blobs_global_window);
            %figure(); imagesc(img_blob_label_global_window); colormap('gray');%---
            
            
            %figure(); imagesc(img_blob_label_global); colormap('gray');%---
            % Removes all blobs inside the window
            for i = 1:blob_count_global_window
                img_global_threshold(find(img_blob_label_global == blobs_global_window(i))) = 0;
                img_blob_label_global(find(img_blob_label_global == blobs_global_window(i))) = 0;
            end
            %figure(); imagesc(img_blob_label_global); colormap('gray');%---
            
            % Count number of blobs in window
            [img_blob_label_window,blob_count_window] = bwlabel(img_temp_window,4);
            %figure(); imagesc(img_blob_label_window); colormap('gray');%---
            
            % Find label value for max intensity point
            blob_label_value = img_blob_label_window(find(img_temp_window == max(img_temp_window(:))));
            
            % Remove all other blobs than max intensity
            img_blob_label_window(find(img_blob_label_window ~= blob_label_value(1))) = 0;
            img_blob_label_window(img_blob_label_window > 0) = 1;
            %img_temp_window(find(img_blob_label_window ~= blob_label_value(1))) = 0;
            
            % Check conditions
            if (MB_window_out_of_bounce == 0)% && (blob_count_global_window <= MB_count_condition) && (blob_count_window <= MB_count_condition)
                
                % Features: MaxIntensity, MeanIntensity, WeigtedCentroid
                MB_blob_features = regionprops(img_blob_label_window, img_temp_window,'MaxIntensity','WeightedCentroid');
                
                % Updata MB
                MB(MB_index).old_pos = MB(MB_index).new_pos;
                MB(MB_index).new_pos = fliplr(round(MB_blob_features.WeightedCentroid)) + [MB_window_coord_search(1,1)-1, MB_window_coord_search(1,2)-1];%[max_y,max_x];
                MB(MB_index).vel = MB(MB_index).new_pos-MB(MB_index).old_pos;
                MB(MB_index).max_int = max_int;
                MB(MB_index).centroid = fliplr(round(MB_blob_features.WeightedCentroid)) + [MB_window_coord_search(1,1)-1, MB_window_coord_search(1,2)-1];
                MB(MB_index).count = [blob_count_global_window blob_count_window];
                MB(MB_index).age = MB(MB_index).age + [0 1 1];

                % Updata MB_log
                MB_log(MB_index).old_pos(MB(MB_index).age(3),:) = MB(MB_index).old_pos;
                MB_log(MB_index).new_pos(MB(MB_index).age(3),:) = MB(MB_index).new_pos;
                MB_log(MB_index).vel(MB(MB_index).age(3),:) = MB(MB_index).vel;
                MB_log(MB_index).max_int(MB(MB_index).age(3)) = MB(MB_index).max_int;
                MB_log(MB_index).centroid(MB(MB_index).age(3),:) = MB(MB_index).centroid;
                MB_log(MB_index).count(MB(MB_index).age(3),:) = MB(MB_index).count;
                MB_log(MB_index).age = MB(MB_index).age;
            else
                % Remove MB from list
                MB(MB_index).state = 0;
                MB_log(MB_index).state = 0;
            end
        else
            % Remove MB from list
            MB(MB_index).state = 0;
            MB_log(MB_index).state = 0;
        end
        % Check next MB from list
    end
    
    while any(img_global_threshold(:))         
        % Find max intensity and its coordinates
        [max_int, max_index] = max(img_global_threshold(:));
        [max_y, max_x] = ind2sub(size(img), max_index);

        % Update window coordinates
        MB_window_coord_search(1,1) = max_y-MB_window_size_search_localization(1);
        MB_window_coord_search(2,1) = max_y+MB_window_size_search_localization(1);
        MB_window_coord_search(1,2) = max_x-MB_window_size_search_localization(2);
        MB_window_coord_search(2,2) = max_x+MB_window_size_search_localization(2);

        % Check for out of bounce
        % y_start
        MB_window_out_of_bounce = 0;
        if MB_window_coord_search(1,1) <= 0
            MB_window_coord_search(1,1) = 1;
            MB_window_out_of_bounce = 1;
        end
        % y_end
        if MB_window_coord_search(2,1) > size(img,1)
            MB_window_coord_search(2,1) = size(img,1);
            MB_window_out_of_bounce = 1;
        end
        % x_start
        if MB_window_coord_search(1,2) <= 0
            MB_window_coord_search(1,2) = 1;
            MB_window_out_of_bounce = 1;
        end
        % x_end
        if MB_window_coord_search(2,2) > size(img,2)
            MB_window_coord_search(2,2) = size(img,2);
            MB_window_out_of_bounce = 1;
        end
        
        % Update clear window coordinates
        MB_window_coord_search_clear(1,1) = max_y-MB_window_size_search_clear(1);
        MB_window_coord_search_clear(2,1) = max_y+MB_window_size_search_clear(1);
        MB_window_coord_search_clear(1,2) = max_x-MB_window_size_search_clear(2);
        MB_window_coord_search_clear(2,2) = max_x+MB_window_size_search_clear(2);
        
        % Check for out of bounce
        % y_start      
        if MB_window_coord_search_clear(1,1) <= 0
            MB_window_coord_search_clear(1,1) = 1;
        end
        % y_end
        if MB_window_coord_search_clear(2,1) > size(img,1)
            MB_window_coord_search_clear(2,1) = size(img,1);
        end
        % x_start
        if MB_window_coord_search_clear(1,2) <= 0
            MB_window_coord_search_clear(1,2) = 1;
        end
        % x_end
        if MB_window_coord_search_clear(2,2) > size(img,2)
            MB_window_coord_search_clear(2,2) = size(img,2);
        end
        
        % Create temp window
        img_temp_window = img(MB_window_coord_search(1,1):MB_window_coord_search(2,1),MB_window_coord_search(1,2):MB_window_coord_search(2,2));
        %figure(); imagesc(img_temp_window); colormap('gray');%---
        
        % Thresholding
        img_temp_window(find(img_temp_window < max_int/MB_window_threshold)) = 0;
        figure(); imagesc(img_temp_window); colormap('gray');%---
        
        % Count number of blobs in window from global img
        img_blob_label_global_window = img_blob_label_global(MB_window_coord_search_clear(1,1):MB_window_coord_search_clear(2,1),MB_window_coord_search_clear(1,2):MB_window_coord_search_clear(2,2));
        blobs_global_window = unique(img_blob_label_global_window);
        if length(blobs_global_window) > 1
            blobs_global_window(1) = [];
        end
        blob_count_global_window = length(blobs_global_window);
        %figure(); imagesc(img_blob_label_global_window); colormap('gray');%---
        
        %figure(); imagesc(img_blob_label_global); colormap('gray');%---
        % Removes all blobs inside the window
        for i = 1:blob_count_global_window
            img_global_threshold(find(img_blob_label_global == blobs_global_window(i))) = 0;
            img_blob_label_global(find(img_blob_label_global == blobs_global_window(i))) = 0;
        end
        %figure(); imagesc(img_blob_label_global); colormap('gray');%---
        
        % Count number of blobs in window
        [img_blob_label_window,blob_count_window] = bwlabel(img_temp_window,4);
        %figure(); imagesc(img_blob_label_window); colormap('gray');%---
        
        % Find label value for max intensity point
        blob_label_value = img_blob_label_window(find(img_temp_window == max(img_temp_window(:))));
        
        % Remove all other blobs than max intensity
        img_blob_label_window(find(img_blob_label_window ~= blob_label_value(1))) = 0;     
        img_blob_label_window(img_blob_label_window > 0) = 1;
        %img_temp_window(find(img_blob_label_window ~= blob_label_value(1))) = 0;     
        
        
        % Check for MB_count condition
        if (MB_window_out_of_bounce == 0)% && (blob_count_global_window <= MB_count_condition) && (blob_count_window <= MB_count_condition)
            % Index update
            MB_index = max([MB_log(:).id]) + 1;
            
            % Features: MaxIntensity, MeanIntensity, WeigtedCentroid
            MB_blob_features = regionprops(img_blob_label_window, img_temp_window,'MaxIntensity','WeightedCentroid');
           
            % Updata MB
            MB(MB_index).state = 1;
            MB(MB_index).old_pos = [0, 0];
            MB(MB_index).new_pos = fliplr(round(MB_blob_features.WeightedCentroid)) + [MB_window_coord_search(1,1)-1, MB_window_coord_search(1,2)-1];%[max_y, max_x];
            MB(MB_index).vel = [0,0];
            MB(MB_index).id = MB_index;
            MB(MB_index).max_int = max_int;
            MB(MB_index).centroid = fliplr(round(MB_blob_features.WeightedCentroid)) + [MB_window_coord_search(1,1)-1, MB_window_coord_search(1,2)-1];
            MB(MB_index).count = [blob_count_global_window blob_count_window];
            MB(MB_index).age = [idx_frame idx_frame 1];

            % Updata MB_log
            MB_log(MB_index).state = 1;
            MB_log(MB_index).old_pos = MB(MB_index).old_pos;
            MB_log(MB_index).new_pos = MB(MB_index).new_pos;
            MB_log(MB_index).vel = MB(MB_index).vel;
            MB_log(MB_index).id = MB_index;
            MB_log(MB_index).max_int = MB(MB_index).max_int;
            MB_log(MB_index).centroid = MB(MB_index).centroid;
            MB_log(MB_index).count = MB(MB_index).count;
            MB_log(MB_index).age = MB(MB_index).age;
        end
    end
end

% Remove first instances of pos and vel due to startup
for MB_index = 1:length(MB)
    MB_log(MB_index).old_pos(1,:) = [];
    MB_log(MB_index).vel(1,:) = [];
end
toc
%%
% 
%%
% img = load_img_contrast(2000,n_bck_grnd,n_bck_grnd_skip,0.5*10^-3);
% [X,Y] = meshgrid(1:size(img,2),1:size(img,1));
% [Xq,Yq] = meshgrid(1:1/19.7:size(img,2),1:1/5.1:size(img,1));
% img = interp2(X,Y,img,Xq,Yq,interpolation_type);
% 
% n_bck_grnd = 10; 
% n_bck_grnd_skip = 100;
% v_MB = 0.5*10^(-3);
% idx_frame_start = 2000;
% for i = idx_frame_start:idx_frame_start+nframe
%     pause(0.3)
%     img = load_img_contrast(i,n_bck_grnd,n_bck_grnd_skip,v_MB);
%     img = interp2(X,Y,img,Xq,Yq,interpolation_type);
%     
% %     % Calculate threshold
%     SEM = std(img(:));
%     ts = tinv(0.99999999,length(img(:))-1);
%     CI = mean(img(:)) + ts*SEM;
%     global_threshold = CI
%     img(img < global_threshold) = 0;
%     
%     figure(1); imagesc(img,[global_threshold,60]); colormap('gray');
%     figure(3); plot(img);
%     ylim([global_threshold 60]);
%     xlim([0 475]);
% end
% 
%%
% outputVideo=VideoWriter('MB_video');
% outputVideo.FrameRate=2;
% open(outputVideo);
% mov(1:50)= struct('cdata',[],'colormap',[]);
%     
for idx_frame = 4026:4046
    
    % Load img
    img = load_img_contrast(idx_frame,n_bck_grnd,n_bck_grnd_skip,v_MB);
    img = interp2(X,Y,img,Xq,Yq,interpolation_type);
    %img = img(15:end,:);
    
    % Calculate threshold
    SEM = std(img(:));               
    ts = tinv(0.9999999999,length(img(:))-1);     
    CI = mean(img(:)) + ts*SEM;
    global_threshold = CI;
    
    img_global_threshold = img;
    img_global_threshold(img_global_threshold < global_threshold) = 0;
    
    %img = load_img_contrast(idx_frame,n_bck_grnd,n_bck_grnd_skip,v_MB);
    pause(0.2)
    figure(4); plot(img_global_threshold);
%     ylim([threshold, threshold + 40]);
%    figure(2); imagesc(img); colormap('gray');
    figure(2); imagesc(img_global_threshold); colormap('gray');    
    xlabel('Lateral [mm]'); ylabel('Axial [mm]'); % title('Micro-Bubble image');
    set(gca,'Xtick',linspace(0,1189,5)); set(gca, 'XTickLabel',linspace(0,12,5));
    set(gca,'Ytick',linspace(0,2489,6)); set(gca, 'YTickLabel',linspace(0,25,6));
    set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
    set(gca, 'PlotBoxAspectRatio',[1 1 1])
%     mov=getframe(gcf);
%     writeVideo(outputVideo,mov.cdata);

    
    
%     norm = max(img_avg(:));
%     limg_avg = 20*log10(img_avg/norm);
%     figure(); imagesc(limg_avg,[-40 0]); colormap(gray);
end
% close(gcf)
% close(outputVideo);

   
%% Find MB from coordinates
temp = [];
coord = [315 486];
for i = 1:size(MB_log,2)
    for j = 1:size(MB_log(i).centroid,1)
        if MB_log(i).centroid(j,1) == coord(2) && MB_log(i).centroid(j,2) == coord(1); 
            temp = [temp i]
            break;
        end
    end
end

%% Plot MB path and frequency response
temp_y = MB_log(temp(1)).centroid(:,1);
temp_x = MB_log(temp(1)).centroid(:,2);
figure(); plot(temp_y)
%figure(); plot(temp_x)te

temp_y_fft = fft(temp_y-mean(temp_y));
figure(); plot(linspace(0,fps/2,round(size(temp_y_fft,1)/2)),abs(temp_y_fft(1:round(size(temp_y_fft,1)/2)))); legend('1','2','3','4','5','6','7');  xlabel('Frequency (Hz)'); ylabel('Magnitude'); % title('Frequency spectrum of axial displacement');
%temp_x_fft = fft(temp_x-mean(temp_x));
%figure(); plot(linspace(0,fps/2,round(size(temp_x_fft,1)/2)),abs(temp_x_fft(1:round(size(temp_x_fft,1)/2)))); legend('1','2','3','4','5','6','7');  xlabel('Frequency (Hz)'); ylabel('Magnitude'); % title('Frequency spectrum of axial displacement');

%%
outputVideo=VideoWriter('contrast_video_clean_log_35');
outputVideo.FrameRate=2;
open(outputVideo);
mov(1:50)= struct('cdata',[],'colormap',[]);
pause(0.3)    
for idx_frame = 4026:4046
    
     % Load img
    img = load_img_contrast(idx_frame,n_bck_grnd,n_bck_grnd_skip,v_MB);
%    img = abs(interp2(X,Y,img,Xq,Yq,interpolation_type));
    %img = img(15:end,:);
    
%     % Calculate threshold
%     SEM = std(img(:));               
%     ts = tinv(0.9999999999,length(img(:))-1);     
%     CI = mean(img(:)) + ts*SEM;
%     global_threshold = CI;
    
%     img_global_threshold = img;
%     img_global_threshold(img_global_threshold < global_threshold) = 0;
    
    %img = load_img_contrast(idx_frame,n_bck_grnd,n_bck_grnd_skip,v_MB);
    pause(0.3)
%     figure(4); plot(img);
%     ylim([0 100]);
    figure(2); imagesc(20*log10(img/max(img(:))), [-35 0]); colormap('gray');
%     figure(2); imagesc(img_global_threshold); colormap('gray');    
    xlabel('Lateral [mm]'); ylabel('Axial [mm]'); % title('Micro-Bubble image');
    set(gca,'Xtick',linspace(0,size(img,2),5)); set(gca, 'XTickLabel',linspace(0,12,5));
    set(gca,'Ytick',linspace(0,size(img,1),6)); set(gca, 'YTickLabel',linspace(0,25,6));
    set(gca, 'DataAspectRatio',[1 4 1]) % set data aspect ratio in zoom box
    set(gca, 'PlotBoxAspectRatio',[1 1 1])
    
    mov=getframe(gcf);
    writeVideo(outputVideo,mov.cdata);

    
    
%     norm = max(img_avg(:));
%     limg_avg = 20*log10(img_avg/norm);
%     figure(); imagesc(limg_avg,[-40 0]); colormap(gray);
end
close(gcf)
close(outputVideo);





