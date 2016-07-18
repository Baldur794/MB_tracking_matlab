%% Conditions
MB_count_condition = 1; % Number of areas within search windows
MB_age_condition = 4; % Tracking length condition
MB_dens_condition_avg = 0; % Density condition average
MB_dens_condition_single = 0; % Density condition single
MB_dens_condition_single_avg = 0; 
MB_window_threshold = 2; % Window Threshold (actual MB_window_threshold = Max_intensity*1/threshold)
MB_window_size_search_localization = [20 20]; % Search window for localization of PSF's
MB_window_size_search_new = [10 10]; % Search window for new MB's
MB_window_size_search_existing = [7 7]; % Search window for ''old'' MB's
MB_window_size_density_avg = [5 3]; % Window for density condition
MB_window_size_density_single = [5 5]; % Window for density condition


% Create annular structuring element for inner marking segmentation
radius_inner = 10;
radius_outer = 30;
center = radius_outer+1;
[W,H] = meshgrid(1:2*radius_outer+1,1:2*radius_outer+1);
strel_ring = (sqrt((W-center).^2 + (H-center).^2) <= radius_outer) & (sqrt((W-center).^2 + (H-center).^2) >= radius_inner); 

% Color matrix for plot
color_plot = [0    0.4470    0.7410
                0.8500    0.3250    0.0980
                0.9290    0.6940    0.1250
                0.4940    0.1840    0.5560
                0.4660    0.6740    0.1880
                0.3010    0.7450    0.9330
                0.6350    0.0780    0.1840];

% img size
img_size = [2418,1183];% flow_behind_BF_res1010
v_MB = 0.5*10^(-3);
threshold = 30;
[X,Y] = meshgrid(1:61,1:475);
[Xq,Yq] = meshgrid(1:1/19.7:61,1:1/5.1:475);

%% Track algorithm
% Start internal timer
tic
MB_window_coord = zeros(2,2); % Coordinates for search window
MB_index = 1; % Index for each MB
MB_window_out_of_bounce = 0; % Checks if search windows is outside image

% initializing MB struct
MB = []; % Current MB
MB_log = []; % MB tracking list
MB(MB_index).state = 0;
MB(MB_index).old_pos = zeros(1,2);
MB(MB_index).new_pos = zeros(1,2);
MB(MB_index).vel = zeros(1,2); 
MB(MB_index).max_int = 0;
MB(MB_index).centroid_I = zeros(1,2);
MB(MB_index).id = 1;
MB(MB_index).count = 0;

idx_seq_start = 4000; % 
idx_seq = idx_seq_start;

img = load_img_contrast(idx_seq_start,0,200,v_MB);
img = img(15:end,:);
img(img < threshold) = 0;
fore_grnd_img = interp2(X,Y,img,Xq,Yq);

% inner maker image
% dilate
fore_grnd_img_inner_marker_dilate = imdilate(fore_grnd_img,strel_ring);

% logical pointwise lowest identifier
pointwise_mask = (fore_grnd_img_inner_marker_dilate > fore_grnd_img);

% pointwise lowest
fore_grnd_img_inner_pointwise_lowest = fore_grnd_img.*pointwise_mask + fore_grnd_img_inner_marker_dilate.*(~pointwise_mask); 

% (original - pointwise lowest)
fore_grnd_img_inner_markers_pointwise_lowest = fore_grnd_img-fore_grnd_img_inner_pointwise_lowest;

% inner markers with original values
fore_grnd_img_inner_markers = fore_grnd_img.* (fore_grnd_img_inner_markers_pointwise_lowest > 0);

% blob labeling of inner markers
[blob_label_inner_marker_image,blob_count_inner_marker] = bwlabel(fore_grnd_img_inner_markers,4);


% Startup:
while any(blob_label_inner_marker_image(:))
    % Find max intensity and its coordinates
    [max_int, max_index] = max(fore_grnd_img_inner_markers(:));
    [max_y, max_x] = ind2sub(size(fore_grnd_img), max_index);
    
    % Update window coordinates
    MB_window_coord(1,1) = max_y-MB_window_size_search_localization(1);
    MB_window_coord(2,1) = max_y+MB_window_size_search_localization(1);
    MB_window_coord(1,2) = max_x-MB_window_size_search_localization(2);
    MB_window_coord(2,2) = max_x+MB_window_size_search_localization(2);
    
    % check for out of bounce
    % y1
    MB_window_out_of_bounce = 0;
    if MB_window_coord(1,1) <= 0
        MB_window_coord(1,1) = 1;
        MB_window_out_of_bounce = 1;
    end
    % y2
    if MB_window_coord(2,1) > size(fore_grnd_img,1)
        MB_window_coord(2,1) = size(fore_grnd_img,1);
        MB_window_out_of_bounce = 1;
    end
    % x1
    if MB_window_coord(1,2) <= 0
        MB_window_coord(1,2) = 1;
        MB_window_out_of_bounce = 1;
    end
    % x2
    if MB_window_coord(2,2) > size(fore_grnd_img,2)
        MB_window_coord(2,2) = size(fore_grnd_img,2);
        MB_window_out_of_bounce = 1;
    end
    
    % Create temp window
    fore_grnd_img_temp_window = fore_grnd_img(MB_window_coord(1,1):MB_window_coord(2,1),MB_window_coord(1,2):MB_window_coord(2,2));
    
    % Thresholding
    fore_grnd_img_temp_window(find(fore_grnd_img_temp_window < max_int/MB_window_threshold)) = 0;

    % Count number of blobs in inner_marker_window
    blob_label_image_inner_marker_window = blob_label_inner_marker_image(MB_window_coord(1,1):MB_window_coord(2,1),MB_window_coord(1,2):MB_window_coord(2,2));
    blobs_in_window = unique(blob_label_image_inner_marker_window);
    blobs_in_window(1) = [];
    blob_count_inner_marker_window = length(blobs_in_window);
    
    % Picks out the blobs inside the blob_label_inner_maker_window
    for i = 1:blob_count_inner_marker_window
        fore_grnd_img_inner_markers(find(blob_label_inner_marker_image == blobs_in_window(i))) = 0;
        blob_label_inner_marker_image(find(blob_label_inner_marker_image == blobs_in_window(i))) = 0;
    end
    
    % Count number of blobs in window
    [blob_label_foregrnd_image_window,blob_count_foregrnd_window] = bwlabel(fore_grnd_img_temp_window,4);
    
    % Find label value for max intensity point
    blob_label_value = blob_label_foregrnd_image_window(find(fore_grnd_img_temp_window == max(fore_grnd_img_temp_window(:))));
    
    % Remove all other blobs than max intensity
    blob_label_foregrnd_image_window(find(blob_label_foregrnd_image_window ~= blob_label_value(1)))  = 0; 
        
    % Check for MB conditions
    if (blob_count_inner_marker_window <= MB_count_condition) && (blob_count_foregrnd_window <= MB_count_condition && MB_window_out_of_bounce == 0)
        % Features I: MaxIntensity, MeanIntensity, WeigtedCentroid
        MB_blob_features_I = regionprops(blob_label_foregrnd_image_window,fore_grnd_img_temp_window,'MaxIntensity','WeightedCentroid');

        % Updata MB
        MB(MB_index).state = 1;
        MB(MB_index).old_pos = [0, 0];
        MB(MB_index).new_pos = [max_y, max_x];
        MB(MB_index).vel = [0,0];
        MB(MB_index).id = MB_index;
        MB(MB_index).max_int = MB_blob_features_I.MaxIntensity;
        MB(MB_index).centroid_I = fliplr(round(MB_blob_features_I.WeightedCentroid)) + [MB_window_coord(1,1)-1, MB_window_coord(1,2)-1];
        MB(MB_index).count = [blob_count_inner_marker_window blob_count_foregrnd_window];
        MB(MB_index).age = [idx_seq idx_seq 1];

        MB_index = MB_index + 1;
    end
end

% Save MB in MB list
MB_log = MB;

% After startup
steps = 1; % number of frames skipped
% Run through all frames
for idx_seq = idx_seq_start + steps:steps:120
    idx_seq
    % Load new frame
    img = load_img_contrast(idx_seq,0,200,v_MB);
    img = img(15:end,:);
    img(img < threshold) = 0;
    fore_grnd_img = interp2(X,Y,img,Xq,Yq);
    
    % Inner maker image
    % Dilate
    fore_grnd_img_inner_marker_dilate = imdilate(fore_grnd_img,strel_ring);
 
    % Logical pointwise lowest identifier
    pointwise_mask = (fore_grnd_img_inner_marker_dilate > fore_grnd_img);
   
    % Pointwise lowest
    fore_grnd_img_inner_pointwise_lowest = fore_grnd_img.*pointwise_mask + fore_grnd_img_inner_marker_dilate.*(~pointwise_mask); 
   
    % (original - pointwise lowest)
    fore_grnd_img_inner_markers_pointwise_lowest = fore_grnd_img-fore_grnd_img_inner_pointwise_lowest;
    
    % Inner markers with original values
    fore_grnd_img_inner_markers = fore_grnd_img.* (fore_grnd_img_inner_markers_pointwise_lowest > 0);
    
    % Blob labeling of inner markers
    [blob_label_inner_marker_image,blob_count_inner_marker] = bwlabel(fore_grnd_img_inner_markers,4);
        
    % Sort MB's after intensity
    MB_sort_index = logical([MB.state]);
    [~,MB_sort_index] = sort([MB.max_int].*MB_sort_index,'descend'); 
    
    % Connect MB's
    for i = 1:sum([MB.state]) 
        % Update MB_index
        MB_index = MB_sort_index(i);
        
        % Update window coordinates
        if MB(MB_index).age(3) == 1
            MB_window_coord(1,1) = MB(MB_index).new_pos(1)-MB_window_size_search_new(1)+MB(MB_index).vel(1);
            MB_window_coord(2,1) = MB(MB_index).new_pos(1)+MB_window_size_search_new(1)+MB(MB_index).vel(1);
            MB_window_coord(1,2) = MB(MB_index).new_pos(2)-MB_window_size_search_new(2)+MB(MB_index).vel(2);
            MB_window_coord(2,2) = MB(MB_index).new_pos(2)+MB_window_size_search_new(2)+MB(MB_index).vel(2);
        else
            MB_window_coord(1,1) = MB(MB_index).new_pos(1)-MB_window_size_search_existing(1)+MB(MB_index).vel(1);
            MB_window_coord(2,1) = MB(MB_index).new_pos(1)+MB_window_size_search_existing(1)+MB(MB_index).vel(1);
            MB_window_coord(1,2) = MB(MB_index).new_pos(2)-MB_window_size_search_existing(2)+MB(MB_index).vel(2);
            MB_window_coord(2,2) = MB(MB_index).new_pos(2)+MB_window_size_search_existing(2)+MB(MB_index).vel(2);
        end
        
        % check for out of bounce
        % y1
        MB_window_out_of_bounce = 0;
        if MB_window_coord(1,1) <= 0
            MB_window_coord(1,1) = 1;
            MB_window_out_of_bounce = 1;
        end
        % y2
        if MB_window_coord(2,1) > size(fore_grnd_img,1)
            MB_window_coord(2,1) = size(fore_grnd_img,1);
            MB_window_out_of_bounce = 1;
        end
        % x1
        if MB_window_coord(1,2) <= 0
            MB_window_coord(1,2) = 1;
            MB_window_out_of_bounce = 1;
        end
        % x2
        if MB_window_coord(2,2) > size(fore_grnd_img,2)
            MB_window_coord(2,2) = size(fore_grnd_img,2);
            MB_window_out_of_bounce = 1;
        end
        
        % Get blob window
        blob_label_image_inner_marker_window = blob_label_inner_marker_image(MB_window_coord(1,1):MB_window_coord(2,1),MB_window_coord(1,2):MB_window_coord(2,2));        

        % Check if any blobs are within window
        if any(blob_label_image_inner_marker_window(:))
            % Count number of blobs in inner_marker_window
            blobs_in_window = unique(blob_label_image_inner_marker_window);
            blobs_in_window(1) = [];
            blob_count_inner_marker_window = length(blobs_in_window);
        
            % Create temp window around the respective MB
            fore_grnd_img_temp_window = fore_grnd_img(MB_window_coord(1,1):MB_window_coord(2,1),MB_window_coord(1,2):MB_window_coord(2,2));

            % Create temp window around the respective MB from inner markers
            fore_grnd_img_inner_markers_window = fore_grnd_img_inner_markers(MB_window_coord(1,1):MB_window_coord(2,1),MB_window_coord(1,2):MB_window_coord(2,2));
            
            % Find max intensity within window
            [max_int, max_index_window] = max(fore_grnd_img_inner_markers_window(:));
            [max_y_window, max_x_window] = ind2sub(size(fore_grnd_img_temp_window), max_index_window);
            max_y = max_y_window + MB_window_coord(1,1)-1;
            max_x = max_x_window + MB_window_coord(1,2)-1;
            
            % Window thresholding
            fore_grnd_img_temp_window(find(fore_grnd_img_temp_window < max_int/MB_window_threshold)) = 0;
            
            % Picks out the blobs inside the blob_label_inner_maker_window
            for i = 1:blob_count_inner_marker_window
                fore_grnd_img_inner_markers(find(blob_label_inner_marker_image == blobs_in_window(i))) = 0;
                blob_label_inner_marker_image(find(blob_label_inner_marker_image == blobs_in_window(i))) = 0;
            end
            
            % Count number of blobs in window
            [blob_label_foregrnd_image_window,blob_count_foregrnd_window] = bwlabel(fore_grnd_img_temp_window,4);
            
            % Find label value for max intensity point
            blob_label_value = blob_label_foregrnd_image_window(find(fore_grnd_img_temp_window == max(fore_grnd_img_temp_window(:))));
            
            % Remove all other blobs than max intensity
            blob_label_foregrnd_image_window(find(blob_label_foregrnd_image_window ~= blob_label_value(1))) = 0; 

            % Check conditions
            if (blob_count_inner_marker_window <= MB_count_condition) && (blob_count_foregrnd_window <= MB_count_condition && MB_window_out_of_bounce == 0)
                % Features I: MaxIntensity, MeanIntensity, WeigtedCentroid
                MB_blob_features_I = regionprops(blob_label_foregrnd_image_window, fore_grnd_img_temp_window,'MaxIntensity','WeightedCentroid');

                % Updata MB
                MB(MB_index).old_pos = MB(MB_index).new_pos;
                MB(MB_index).new_pos = [max_y,max_x];
                MB(MB_index).vel = MB(MB_index).new_pos-MB(MB_index).old_pos;
                MB(MB_index).max_int = MB_blob_features_I.MaxIntensity;
                MB(MB_index).centroid_I = fliplr(round(MB_blob_features_I.WeightedCentroid)) + [MB_window_coord(1,1)-1, MB_window_coord(1,2)-1];
                MB(MB_index).count = [blob_count_inner_marker_window blob_count_foregrnd_window];
                MB(MB_index).age = MB(MB_index).age + [0 steps 1];

                % Updata MB_log
                MB_log(MB_index).old_pos(MB(MB_index).age(3),:) = MB(MB_index).old_pos;
                MB_log(MB_index).new_pos(MB(MB_index).age(3),:) = MB(MB_index).new_pos;
                MB_log(MB_index).vel(MB(MB_index).age(3),:) = MB(MB_index).vel;
                MB_log(MB_index).max_int(MB(MB_index).age(3)) = MB(MB_index).max_int;
                MB_log(MB_index).centroid_I(MB(MB_index).age(3),:) = MB(MB_index).centroid_I;
                MB_log(MB_index).count(MB(MB_index).age(3),:) = MB(MB_index).count;
                MB_log(MB_index).age = MB(MB_index).age;
            else
                % Remove MB from list
                MB_log(MB_index).state = 0;
            end
        else
            % Remove MB from list
            MB_log(MB_index).state = 0;
        end
        % Check next MB from list
    end
    
    % Seach for new MB's
    while any(blob_label_inner_marker_image(:))         
        % Find max intensity and its coordinates
        [max_int, max_index] = max(fore_grnd_img_inner_markers(:));
        [max_y, max_x] = ind2sub(size(fore_grnd_img), max_index);

        % Update window coordinates
        MB_window_coord(1,1) = max_y-MB_window_size_search_localization(1);
        MB_window_coord(2,1) = max_y+MB_window_size_search_localization(1);
        MB_window_coord(1,2) = max_x-MB_window_size_search_localization(2);
        MB_window_coord(2,2) = max_x+MB_window_size_search_localization(2);

        % Check for out of bounce
        % y1
        MB_window_out_of_bounce = 0;
        if MB_window_coord(1,1) <= 0
            MB_window_coord(1,1) = 1;
            MB_window_out_of_bounce = 1;
        end
        % y2
        if MB_window_coord(2,1) > size(fore_grnd_img,1)
            MB_window_coord(2,1) = size(fore_grnd_img,1);
            MB_window_out_of_bounce = 1;
        end
        % x1
        if MB_window_coord(1,2) <= 0
            MB_window_coord(1,2) = 1;
            MB_window_out_of_bounce = 1;
        end
        % x2
        if MB_window_coord(2,2) > size(fore_grnd_img,2)
            MB_window_coord(2,2) = size(fore_grnd_img,2);
            MB_window_out_of_bounce = 1;
        end

        % Create temp window
        fore_grnd_img_temp_window = fore_grnd_img(MB_window_coord(1,1):MB_window_coord(2,1),MB_window_coord(1,2):MB_window_coord(2,2));

        % Window thresholding
        fore_grnd_img_temp_window(find(fore_grnd_img_temp_window < max_int/MB_window_threshold)) = 0;

        % Count number of blobs in inner_marker_window
        blob_label_image_inner_marker_window = blob_label_inner_marker_image(MB_window_coord(1,1):MB_window_coord(2,1),MB_window_coord(1,2):MB_window_coord(2,2));
        blobs_in_window = unique(blob_label_image_inner_marker_window);
        blobs_in_window(1) = [];
        blob_count_inner_marker_window = length(blobs_in_window);

        % Picks out the blobs inside the blob_label_inner_maker_window
        for i = 1:blob_count_inner_marker_window
            fore_grnd_img_inner_markers(find(blob_label_inner_marker_image == blobs_in_window(i))) = 0;
            blob_label_inner_marker_image(find(blob_label_inner_marker_image == blobs_in_window(i))) = 0;
        end

        % Count number of blobs in window
        [blob_label_foregrnd_image_window,blob_count_foregrnd_window] = bwlabel(fore_grnd_img_temp_window,4);
        
        % Find label value for max intensity point
        blob_label_value = blob_label_foregrnd_image_window(find(fore_grnd_img_temp_window == max(fore_grnd_img_temp_window(:))));
        
        % Remove all other blobs than max intensity
        blob_label_foregrnd_image_window(find(blob_label_foregrnd_image_window ~= blob_label_value(1)))  = 0; 
        
        % Check for MB_count condition
        if (blob_count_inner_marker_window <= MB_count_condition) && (blob_count_foregrnd_window <= MB_count_condition && MB_window_out_of_bounce == 0)
            % Index update
            MB_index = max([MB_log(:).id]) + 1;
            
            % Picks out the blob with max intensity
            blob_label_image_extract_blobs = zeros(size(blob_label_inner_marker_image));
            blob_label_image_extract_blobs(find(blob_label_inner_marker_image == blob_label_inner_marker_image(max_y,max_x))) = 1; % max_y, max_x legit here..

            % Features I: MaxIntensity, MeanIntensity, WeigtedCentroid
            MB_blob_features_I = regionprops(blob_label_foregrnd_image_window,fore_grnd_img_temp_window,'MaxIntensity','WeightedCentroid');

            % Set respective blobs to zero
            fore_grnd_img = fore_grnd_img.*(~blob_label_image_extract_blobs);

            % Updata MB
            MB(MB_index).state = 1;
            MB(MB_index).old_pos = [0, 0];
            MB(MB_index).new_pos = [max_y, max_x];
            MB(MB_index).vel = [0,0];
            MB(MB_index).id = MB_index;
            MB(MB_index).max_int = MB_blob_features_I.MaxIntensity;
            MB(MB_index).centroid_I = fliplr(round(MB_blob_features_I.WeightedCentroid)) + [MB_window_coord(1,1)-1, MB_window_coord(1,2)-1];
            MB(MB_index).count = [blob_count_inner_marker_window blob_count_foregrnd_window];
            MB(MB_index).age = [idx_seq idx_seq 1];

            % Updata MB_log
            MB_log(MB_index).state = 1;
            MB_log(MB_index).old_pos = MB(MB_index).old_pos;
            MB_log(MB_index).new_pos = MB(MB_index).new_pos;
            MB_log(MB_index).vel = MB(MB_index).vel;
            MB_log(MB_index).id = MB_index;
            MB_log(MB_index).max_int = MB(MB_index).max_int;
            MB_log(MB_index).centroid_I = MB(MB_index).centroid_I;
            MB_log(MB_index).count = MB(MB_index).count;
            MB_log(MB_index).age = MB(MB_index).age;
        end
    end
end

% Running time
total_time = toc % total time

% Remove first instances of pos and vel due to startup
for MB_index = 1:length(MB)
    MB_log(MB_index).old_pos(1,:) = [];
    MB_log(MB_index).vel(1,:) = [];
end