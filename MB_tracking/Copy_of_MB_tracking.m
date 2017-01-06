log_compression = @(x) 20*log10(x/max(x(:)));

%% Tracking parameters
% Tracking parameters
MB_window_coord_search = zeros(2,2); % Coordinates for search window
MB_window_coord_search_clear = zeros(2,2); % Coordinates for search clear window

MB_window_size_search_localization = [30 30]; % [y,x] Search window for localization of PSF's
MB_window_size_search_new = [3 3]; % [y,x] Search window for new MB's
MB_window_size_search_existing = [3 3]; % [y,x] Search window for ''old'' MB's
MB_window_size_search_clear = [40 40]; % [y,x] Search window for ''old'' MB's
MB_window_threshold = 1.3; % Window Threshold (actual MB_window_threshold = Max_intensity*1/threshold)


idx_frame_start = 500;
nframe = 100;

% initializing MB struct

MB = []; % Current MB
MB.state = 0;
MB.old_pos = zeros(1,2);
MB.new_pos = zeros(1,2);
MB.vel = zeros(1,2); 
MB.max_int = 0;
MB.centroid = zeros(1,2);
MB.area = 0;
MB.eccentricity = 0;
MB.orientation = 0;
MB.perimeter = 0;
MB.id = 0;
MB.count = 0;
MB_log = MB;

MB_index = 0; % Index for each MB
MB_window_out_of_bounce = 0; % Checks if search windows is outside image

tic
for idx_frame=idx_frame_start:idx_frame_start+nframe
    idx_frame
    
    % Load img
    [PI_img, time_stamp, line_count] = load_img(idx_frame, img_type, folderName, img_resolution, area_of_interest);;

    % Calculate threshold
    SEM = std(PI_img(:));
    ts = tinv(0.9999,length(PI_img(:))-1);     
    CI = mean(PI_img(:)) + ts*SEM;
    global_threshold = CI;

    
    % Threshold
    img_global_threshold = PI_img;
    img_global_threshold(img_global_threshold < global_threshold) = 0;
%     figure(); imagesc(img_global_threshold); colormap('gray');%---
    
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
        if MB(MB_index).age(3) == 1
            MB_window_coord_search(1,1) = MB(MB_index).new_pos(1)-MB_window_size_search_new(1);
            MB_window_coord_search(2,1) = MB(MB_index).new_pos(1)+MB_window_size_search_new(1);
            MB_window_coord_search(1,2) = MB(MB_index).new_pos(2)-MB_window_size_search_new(2);
            MB_window_coord_search(2,2) = MB(MB_index).new_pos(2)+MB_window_size_search_new(2);
        else
            MB_window_coord_search(1,1) = MB(MB_index).new_pos(1)-MB_window_size_search_existing(1);%+MB(MB_index).vel(1);
            MB_window_coord_search(2,1) = MB(MB_index).new_pos(1)+MB_window_size_search_existing(1);%+MB(MB_index).vel(1);
            MB_window_coord_search(1,2) = MB(MB_index).new_pos(2)-MB_window_size_search_existing(2);%+MB(MB_index).vel(2);
            MB_window_coord_search(2,2) = MB(MB_index).new_pos(2)+MB_window_size_search_existing(2);%+MB(MB_index).vel(2);
        end
        
        % Check for out of bounce
        % y_start
        MB_window_out_of_bounce = 0;
        if MB_window_coord_search(1,1) <= 0
            MB_window_coord_search(1,1) = 1;
            MB_window_out_of_bounce = 1;
        end
        % y_end
        if MB_window_coord_search(2,1) > size(PI_img,1)
            MB_window_coord_search(2,1) = size(PI_img,1);
            MB_window_out_of_bounce = 1;
        end
        % x_start
        if MB_window_coord_search(1,2) <= 0
            MB_window_coord_search(1,2) = 1;
            MB_window_out_of_bounce = 1;
        end
        % x_end
        if MB_window_coord_search(2,2) > size(PI_img,2)
            MB_window_coord_search(2,2) = size(PI_img,2);
            MB_window_out_of_bounce = 1;
        end       
        
        % Create temp window
        img_temp_window = img_global_threshold(MB_window_coord_search(1,1):MB_window_coord_search(2,1),MB_window_coord_search(1,2):MB_window_coord_search(2,2));
        
         % Check if any blobs are within window
        if any(img_temp_window(:))                            
            % Find max intensity within window
            [max_int, max_index] = max(img_temp_window(:));
            [max_y_window, max_x_window] = ind2sub(size(img_temp_window), max_index);
            max_y = max_y_window + MB_window_coord_search(1,1)-1;
            max_x = max_x_window + MB_window_coord_search(1,2)-1;
                                            
            % Update window coordinates
            MB_window_coord_search(1,1) = max_y-MB_window_size_search_localization(1);
            MB_window_coord_search(2,1) = max_y+MB_window_size_search_localization(1);
            MB_window_coord_search(1,2) = max_x-MB_window_size_search_localization(2);
            MB_window_coord_search(2,2) = max_x+MB_window_size_search_localization(2);
            
            % Update clear window coordinates
            MB_window_coord_search_clear(1,1) = max_y-MB_window_size_search_clear(1);
            MB_window_coord_search_clear(2,1) = max_y+MB_window_size_search_clear(1);
            MB_window_coord_search_clear(1,2) = max_x-MB_window_size_search_clear(2);
            MB_window_coord_search_clear(2,2) = max_x+MB_window_size_search_clear(2);
                        
            % Check for out of bounce
            % No need to check for MB_window_coord_search if MB_window_coord_search < MB_window_coord_search_clear... 
            % y_start
            MB_window_out_of_bounce = 0;
            if MB_window_coord_search(1,1) <= 0 || MB_window_coord_search_clear(1,1) <= 0
                MB_window_coord_search(1,1) = 1;
                MB_window_coord_search_clear(1,1) = 1;
                MB_window_out_of_bounce = 1;
            end
            % y_end
            if MB_window_coord_search(2,1) > size(PI_img,1) || MB_window_coord_search_clear(2,1) > size(PI_img,1)
                MB_window_coord_search(2,1) = size(PI_img,1);
                MB_window_coord_search_clear(2,1) = size(PI_img,1);
                MB_window_out_of_bounce = 1;
            end
            % x_start
            if MB_window_coord_search(1,2) <= 0 || MB_window_coord_search_clear(1,2) <= 0
                MB_window_coord_search(1,2) = 1;
                MB_window_coord_search_clear(1,2) = 1;
                MB_window_out_of_bounce = 1;
            end
            % x_end
            if MB_window_coord_search(2,2) > size(PI_img,2) || MB_window_coord_search_clear(2,2) > size(PI_img,2)
                MB_window_coord_search(2,2) = size(PI_img,2);
                MB_window_coord_search_clear(2,2) = size(PI_img,2);
                MB_window_out_of_bounce = 1;
            end
                      
            % Create temp window
            img_temp_window = PI_img(MB_window_coord_search(1,1):MB_window_coord_search(2,1),MB_window_coord_search(1,2):MB_window_coord_search(2,2));
%             figure(); imagesc(img_temp_window); colormap('gray');%---
            
            % Local thresholding
            img_temp_window(find(img_temp_window < max_int/MB_window_threshold | img_temp_window > max_int)) = 0;
%             figure(); imagesc(img_temp_window); colormap('gray');%---
            
            % Count number of blobs in window from global img
            img_blob_label_global_window = img_blob_label_global(MB_window_coord_search_clear(1,1):MB_window_coord_search_clear(2,1),MB_window_coord_search_clear(1,2):MB_window_coord_search_clear(2,2));
            blobs_global_window = unique(img_blob_label_global_window);
            blobs_global_window(1) = [];
            blob_count_global_window = length(blobs_global_window);
%             figure(); imagesc(img_blob_label_global_window); colormap('gray');%---
            
            
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
            blob_label_value = img_blob_label_window(find(img_temp_window == max_int));%max(img_temp_window(:))));
            
            % Remove all other blobs than max intensity
            img_blob_label_window(find(img_blob_label_window ~= blob_label_value(1))) = 0;
            img_blob_label_window(img_blob_label_window > 0) = 1;
            %img_temp_window(find(img_blob_label_window ~= blob_label_value(1))) = 0;
            
            % Check conditions
            if (MB_window_out_of_bounce == 0)% && (blob_count_global_window <= MB_count_condition) && (blob_count_window <= MB_count_condition)
                
                % Features: MaxIntensity, MeanIntensity, WeigtedCentroid
                MB_blob_features = regionprops(img_blob_label_window, img_temp_window,'MaxIntensity','WeightedCentroid');
                
                % Features: Area, Eccentricity, Orientation, Perimeter
                MB_blob_features_bw = regionprops(img_blob_label_window,'Area','Eccentricity','Orientation','Perimeter');     
                
                % Updata MB
                MB(MB_index).old_pos = MB(MB_index).new_pos;
                MB(MB_index).new_pos = fliplr(round(MB_blob_features.WeightedCentroid)) + [MB_window_coord_search(1,1)-1, MB_window_coord_search(1,2)-1];%[max_y,max_x];
                MB(MB_index).vel = MB(MB_index).new_pos-MB(MB_index).old_pos;
                MB(MB_index).max_int = max_int;
                MB(MB_index).centroid = fliplr(round(MB_blob_features.WeightedCentroid)) + [MB_window_coord_search(1,1)-1, MB_window_coord_search(1,2)-1];
                MB(MB_index).area = MB_blob_features_bw.Area;
                MB(MB_index).eccentricity = MB_blob_features_bw.Eccentricity;
                MB(MB_index).orientation = MB_blob_features_bw.Orientation;
                MB(MB_index).perimeter = MB_blob_features_bw.Perimeter;
                MB(MB_index).count = [blob_count_global_window blob_count_window];
                MB(MB_index).age = MB(MB_index).age + [0 1 1];

                % Updata MB_log
                MB_log(MB_index).old_pos(MB(MB_index).age(3),:) = MB(MB_index).old_pos;
                MB_log(MB_index).new_pos(MB(MB_index).age(3),:) = MB(MB_index).new_pos;
                MB_log(MB_index).vel(MB(MB_index).age(3),:) = MB(MB_index).vel;
                MB_log(MB_index).max_int(MB(MB_index).age(3)) = MB(MB_index).max_int;
                MB_log(MB_index).centroid(MB(MB_index).age(3),:) = MB(MB_index).centroid;               
                MB_log(MB_index).area(MB(MB_index).age(3)) = MB(MB_index).area;
                MB_log(MB_index).eccentricity(MB(MB_index).age(3)) = MB(MB_index).eccentricity;
                MB_log(MB_index).orientation(MB(MB_index).age(3)) = MB(MB_index).orientation;
                MB_log(MB_index).perimeter(MB(MB_index).age(3)) = MB(MB_index).perimeter;              
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
        [max_y, max_x] = ind2sub(size(PI_img), max_index);

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
        if MB_window_coord_search(2,1) > size(PI_img,1)
            MB_window_coord_search(2,1) = size(PI_img,1);
            MB_window_out_of_bounce = 1;
        end
        % x_start
        if MB_window_coord_search(1,2) <= 0
            MB_window_coord_search(1,2) = 1;
            MB_window_out_of_bounce = 1;
        end
        % x_end
        if MB_window_coord_search(2,2) > size(PI_img,2)
            MB_window_coord_search(2,2) = size(PI_img,2);
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
        if MB_window_coord_search_clear(2,1) > size(PI_img,1)
            MB_window_coord_search_clear(2,1) = size(PI_img,1);
        end
        % x_start
        if MB_window_coord_search_clear(1,2) <= 0
            MB_window_coord_search_clear(1,2) = 1;
        end
        % x_end
        if MB_window_coord_search_clear(2,2) > size(PI_img,2)
            MB_window_coord_search_clear(2,2) = size(PI_img,2);
        end
        
        % Create temp window
        img_temp_window = PI_img(MB_window_coord_search(1,1):MB_window_coord_search(2,1),MB_window_coord_search(1,2):MB_window_coord_search(2,2));
%         figure(); imagesc(img_temp_window); colormap('gray');%---
        
        % Local thresholding
        img_temp_window(find(img_temp_window < max_int/MB_window_threshold | img_temp_window > max_int)) = 0;
%         figure(); imagesc(img_temp_window); colormap('gray');%---
%         
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
            
            % Features: Area, Eccentricity, Orientation, Perimeter
            MB_blob_features_bw = regionprops(img_blob_label_window,'Area','Eccentricity','Orientation','Perimeter');
            
            % Updata MB
            MB(MB_index).state = 1;
            MB(MB_index).old_pos = [0, 0];
            MB(MB_index).new_pos = fliplr(round(MB_blob_features.WeightedCentroid)) + [MB_window_coord_search(1,1)-1, MB_window_coord_search(1,2)-1];%[max_y, max_x];
            MB(MB_index).vel = [0,0];
            MB(MB_index).id = MB_index;
            MB(MB_index).max_int = max_int;
            MB(MB_index).centroid = fliplr(round(MB_blob_features.WeightedCentroid)) + [MB_window_coord_search(1,1)-1, MB_window_coord_search(1,2)-1];
            MB(MB_index).area = MB_blob_features_bw.Area;
            MB(MB_index).eccentricity = MB_blob_features_bw.Eccentricity;
            MB(MB_index).orientation = MB_blob_features_bw.Orientation;
            MB(MB_index).perimeter = MB_blob_features_bw.Perimeter;
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
            MB_log(MB_index).area = MB(MB_index).area;
            MB_log(MB_index).eccentricity = MB(MB_index).eccentricity;
            MB_log(MB_index).orientation = MB(MB_index).orientation;
            MB_log(MB_index).perimeter = MB(MB_index).perimeter;           
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

%% Find MB from coordinates
MB_index_valid = [];
coord = [365 989];
for i = 1:size(MB_index_filter_age,2)
    for j = 1:size(MB_log(MB_index_filter_age(i)).centroid,1)
        if MB_log(MB_index_filter_age(i)).centroid(j,1) == coord(2) && MB_log(MB_index_filter_age(i)).centroid(j,2) == coord(1); 
            MB_index_valid = [MB_index_valid MB_index_filter_age(i)]
            MB_log(MB_index_filter_age(i)).age
            break
        end
    end
end


%% Plot MB path and frequency response
for i = 1:size(MB_index_valid,2)
    temp_y = MB_log(MB_index_valid(i)).centroid(:,1);
    temp_x = MB_log(MB_index_valid(i)).centroid(:,2);
    figure(); plot(temp_y);
    figure(); plot(temp_x)
    temp_int = MB_log(MB_index_valid(i)).max_int(:);
    temp_int = sum(MB_log(MB_index_valid(i)).max_int(:))/MB_log(MB_index_valid(i)).age(3)
    %temp_y_fft = fft(temp_y-mean(temp_y));
    %figure(); plot(abs(temp_y_fft));
    % %temp_x_fft = fft(temp_x-mean(temp_x));
    %figure(); plot(linspace(0,fps/2,round(size(temp_x_fft,1)/2)),abs(temp_x_fft(1:round(size(temp_x_fft,1)/2)))); legend('1','2','3','4','5','6','7');  xlabel('Frequency (Hz)'); ylabel('Magnitude'); % title('Frequency spectrum of axial displacement');
end
%% Axial
MB_index_valid = MB_index_filter_dir;
age_max = 0;
age_min = MB_log(MB_index_valid(1)).age(1);
for i = 1:size(MB_index_valid,2)
    MB_log(MB_index_valid(i)).age;
    if MB_log(MB_index_valid(i)).age(2) > age_max
        age_max = MB_log(MB_index_valid(i)).age(2);
    end
    if MB_log(MB_index_valid(i)).age(1) < age_min
        age_min = MB_log(MB_index_valid(i)).age(1);
    end
end
age_min = idx_comp_sync; %--------------------


% Axial
axial_mov = zeros(size(MB_index_valid,2),age_max);
temp_age = zeros(1,age_max);
for i = 1:size(MB_index_valid,2)
    temp = (MB_log(MB_index_valid(i)).centroid(:,1))-(MB_log(MB_index_valid(i)).centroid(end,1)-(MB_log(MB_index_valid(i)).centroid(1,1)))/MB_log(MB_index_valid(i)).age(3)*(1:MB_log(MB_index_valid(i)).age(3))';
    temp = temp-mean(temp);
    axial_mov(i,MB_log(MB_index_valid(i)).age(1):MB_log(MB_index_valid(i)).age(2)) = temp/max(abs(temp));
    temp_age(1,MB_log(MB_index_valid(i)).age(1):MB_log(MB_index_valid(i)).age(2)) = temp_age(1,MB_log(MB_index_valid(i)).age(1):MB_log(MB_index_valid(i)).age(2)) + 1;
end
temp_age(temp_age == 0) = 1;

axial_mov = axial_mov(:,age_min:end);
temp_age = temp_age(1,age_min:end);
axial_mov = sum(axial_mov)./temp_age;
% figure(); plot(axial_mov');
axial_mov_fft = fft(axial_mov);
figure(); plot(abs(axial_mov_fft));
%%
% Lateral
n = 20;
b_h = fir1(n,0.03,'high');

lateral_mov = zeros(size(MB_index_valid,2),age_max);
temp_age = zeros(1,age_max);
for i = 1:size(MB_index_valid,2)
    temp = (MB_log(MB_index_valid(i)).centroid(:,2))-(MB_log(MB_index_valid(i)).centroid(end,2)-(MB_log(MB_index_valid(i)).centroid(1,2)))/MB_log(MB_index_valid(i)).age(3)*(1:MB_log(MB_index_valid(i)).age(3))';
    temp = temp-mean(temp);
    temp = filter(b_h,1,temp,[],1);
    temp(1:n/2-1) = 0;
    temp = circshift(temp,-n/2+1);
    lateral_mov(i,MB_log(MB_index_valid(i)).age(1):MB_log(MB_index_valid(i)).age(2)) = temp/max(abs(temp));
    temp_age(1,MB_log(MB_index_valid(i)).age(1):MB_log(MB_index_valid(i)).age(2)) = temp_age(1,MB_log(MB_index_valid(i)).age(1):MB_log(MB_index_valid(i)).age(2)) + 1;
end
temp_age(temp_age == 0) = 1;

lateral_mov = lateral_mov(:,age_min:end);
temp_age = temp_age(1,age_min:end);
lateral_mov = sum(lateral_mov)./temp_age;
% figure(); plot(lateral_mov');
lateral_mov_fft = fft(lateral_mov);
figure(); plot(abs(lateral_mov_fft));



%% 


figure(2);
% subplot(1,2,1)
PI_img = load_img(518);
imagesc(lateral_axis*1000,depth_axis*1000,log_compression(abs(PI_img)),[-30 0])
xlabel('lateral [mm]')
ylabel('depth [mm]')
axis image
colormap(gray(255))

subplot(1,2,2)
PI_img = load_img(518);
imagesc(lateral_axis*1000,depth_axis*1000,log_compression(abs(PI_img)),[-30 0])
xlabel('lateral [mm]')
ylabel('depth [mm]')
axis image
colormap(gray(255))

%% Make video
% Set image axis:
figure(3); 
PI_img = load_img(12000);


outputVideo=VideoWriter('test5_vid');
outputVideo.FrameRate=3;
open(outputVideo);
mov(1:50)= struct('cdata',[],'colormap',[]);
pause(0.3)
frame_timestamp_linecount = [];
idx = 0;
for idx_frame = 500:700
    
    idx = idx + 1;
    [PI_img, time_stamp, line_count] = load_img(idx_frame);
    frame_timestamp_linecount(idx,:) = [time_stamp, line_count, idx_frame];
    if idx>1
        1/(diff(frame_timestamp_linecount(idx-1:idx,1),1,1)*1e-6)
    end
    
    figure(3);
    imagesc(lateral_axis*1000,depth_axis*1000,log_compression(abs(PI_img)),[-30 0])
%      imagesc(lateral_axis*1000,depth_axis*1000,abs(PI_img),[100 400])
%      imagesc(abs(PI_img),[100 400])
    
    xlabel('lateral [mm]')
    ylabel('depth [mm]')
    axis image
    colormap(gray(255))
    pause(0.1)
    mov=getframe(gcf);
    writeVideo(outputVideo,mov.cdata);
end
close(gcf)
close(outputVideo);

framerate = 1./(diff(frame_timestamp_linecount(:,1),1,1)*1e-6);
figure(2);plot(framerate)






