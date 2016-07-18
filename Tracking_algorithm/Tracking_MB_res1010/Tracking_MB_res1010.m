%% Conditions
MB_count_condition = 1; % 10
MB_age_condition = 1;
MB_dens_condition_avg = 0;
MB_dens_condition_single = 0;
MB_window_size_search_new = [20 20]; % half search window size
MB_window_size_search_existing = [10 10]; % half search window size

% Threshold
threshold = 1.5;

% Create annular structuring element for inner marking segmentation
radius_inner = 4;
radius_outer = 20;
center = radius_outer+1;
[W,H] = meshgrid(1:2*radius_outer+1,1:2*radius_outer+1);
strel_ring = (sqrt((W-center).^2 + (H-center).^2) <= radius_outer) & (sqrt((W-center).^2 + (H-center).^2) >= radius_inner); 

% Create structuring element for erosion of small individual pixes
strel_small_ring = strel('disk',1);

% Color matrix for plot
color_plot = [0    0.4470    0.7410
                0.8500    0.3250    0.0980
                0.9290    0.6940    0.1250
                0.4940    0.1840    0.5560
                0.4660    0.6740    0.1880
                0.3010    0.7450    0.9330
                0.6350    0.0780    0.1840];

% Folderpath
addpath('../'); % for show_img.m and load_data.m

% img size
img_size = [999,301];% flow_behind_BF_res1010
%img_size = [4250,401];% flow_behind_BF_res1010_full_new

% Known parameters  %%

% Axial parameters %
start_depth= 0.028; % Starting depth to beamform [meters]
dr_depth = 10*10^(-6); % Original resolution on the depth axis
end_depth= start_depth + size(img_size,1)*dr_depth;   % Ending depth to beamform [meters]

% Lateral parameters %
start_lateral = -0.009;  % Start with leftmost element
dr_lateral = 10*10^(-6);%1.488281250000001e-04; % Original resolution on the lateral axis
end_lateral = start_lateral + size(img_size,2)*dr_lateral; % End with rightmost element

% Frame rate seq -> seq
frames_rate = 1000/6; % Emission freq = 1kHz




%%



% %% Script to gain info about MB by visuel inspection.
% MB_array = []; % contains all chosen squares
% MB_size_array = []; % list of all chosen square sizes
% MB_list = [] % list of all MB values 
% 
% k = 0;
% idx_seq = 60
% while k == 0 % repeat until keyboard click
%     % load img
%     fore_grnd_img = load_img_tracking_MB_res1010(idx_seq);
%     
%     % set global threshold
%     threshold = max(fore_grnd_img(:))/2;
%     
%     % display the respective img
%     figure, imagesc(fore_grnd_img,[threshold threshold*2]), colormap(gray); 
%     
%     % makes zoom available
%     zoom on;
%     pause()
%     zoom off; 
%     
%     % top left/buttom right coordinates of chosen square
%     MB_square_coord = round(ginput(2));
%    
%     % extract the square within the chosen coordinates
%     MB_square = fore_grnd_img(MB_square_coord(1,2):MB_square_coord(2,2),MB_square_coord(1,1):MB_square_coord(2,1));
%     
%     % log square size
%     MB_size_array = [MB_size_array;[size(MB_square,1),size(MB_square,2)]];
%     
%     % 3d array with all chosen square
%     MB_array(1:size(MB_square,1),1:size(MB_square,2),size(MB_array,3)+1) = MB_square;
%     
%     % list of all MB values 
%     MB_list(size(MB_list,1)+1:size(MB_list,1)+size(MB_square,1)*size(MB_square,2)) = MB_square(:);  
%     
%     k = waitforbuttonpress; % mouse click -> 0, keyboard click -> 1
%     close
%     idx_seq = idx_seq + 5;
% end
% MB_array(:,:,1) = [];
% 
% %% New Script to manually track by visuel inspection.
% MB_index = 1;
% 
% % initializing MB struct
% MB = [];
% MB(MB_index).old_pos = zeros(1,2);
% MB(MB_index).new_pos = zeros(1,2);
% MB(MB_index).vel = zeros(1,2);
% MB(MB_index).age = [1 1 1];
% 
% MB_log = MB;
% 
% k = 0;
% idx_seq = 60;
% while k == 0 % repeat until keyboard click
%     % load img
%     fore_grnd_img = load_img_tracking_MB_res1010(idx_seq);
%     
%     % set global threshold
%     threshold = max(fore_grnd_img(:))/2;
%     
%     %figure, imagesc(abs_data_cropped(:,:,idx_seq)), colormap(gray); % display the respective img
%     figure, imagesc(fore_grnd_img,[threshold threshold*2]), colormap(gray);
%     zoom on; % makes zoom available
%     pause()
%     zoom off;    
%     
%     MB_square_center = round(ginput(1)); % contains top left/buttom right coordinates of chosen square
%     MB(MB_index).old_pos = MB(MB_index).new_pos;
%     MB(MB_index).new_pos = [MB_square_center(2), MB_square_center(1)];
%     MB(MB_index).vel = MB(MB_index).new_pos-MB(MB_index).old_pos;
%     MB(MB_index).age = MB(MB_index).age + [0 1 1];
%     
%     MB_log(MB_index).old_pos(MB(MB_index).age(3),:) = MB(MB_index).old_pos;
%     MB_log(MB_index).new_pos(MB(MB_index).age(3),:) = MB(MB_index).new_pos;
%     MB_log(MB_index).vel(MB(MB_index).age(3),:) = MB(MB_index).vel;
%     MB_log(MB_index).age = MB(MB_index).age;
%     
%     k = waitforbuttonpress; % mouse click -> 0, keyboard click -> 1
% 
%     close
%     idx_seq = idx_seq + 3
% end
% 
% MB_log(MB_index).old_pos(1:2,:) = [];
% MB_log(MB_index).vel(1:2,:) = [];

%% Track algorithm
% Start internal timer
tic
MB_window_coord = zeros(2,2); % coordinates of search window
MB_index = 1; % index for each MB

% initializing MB struct
MB = [];
MB_log = [];
MB(MB_index).state = 0;
MB(MB_index).old_pos = zeros(1,2);
MB(MB_index).new_pos = zeros(1,2);
MB(MB_index).vel = zeros(1,2);
MB(MB_index).max_int = 0;
MB(MB_index).mean_int = 0;
MB(MB_index).centroid_I = zeros(1,2);
MB(MB_index).centroid_BW = zeros(1,2);
MB(MB_index).age = zeros(1,3);
MB(MB_index).id = 1;
MB(MB_index).area = 0;
MB(MB_index).bounding_box = zeros(1,4);
MB(MB_index).eccentricity = 0;
MB(MB_index).extent = 0;
MB(MB_index).major_ax = 0;
MB(MB_index).minor_ax = 0;
MB(MB_index).orientation = 0; 
MB(MB_index).perimeter = 0;
MB(MB_index).count = 0;

idx_seq_start = 101; % frame index-----------------------------------------
idx_seq = idx_seq_start;
fore_grnd_img = load_img_tracking_MB_res1010(idx_seq); % creates a temporary img
%fore_grnd_img = interp2(fore_grnd_img,Xq,Yq); % Interpolates with the factor interp2_x and interp2_y

% inner maker image
% dilate
fore_grnd_img_inner_marker_dilate = imdilate(fore_grnd_img,strel_ring);

% logical pointwise lowest identifier
pointwise_mask = (fore_grnd_img_inner_marker_dilate > fore_grnd_img);

% pointwise lowest
fore_grnd_img_inner_pointwise_lowest = fore_grnd_img.*pointwise_mask + fore_grnd_img_inner_marker_dilate.*(~pointwise_mask); 

% (original - pointwise lowest)
fore_grnd_img_inner_markers_pointwise_lowest = fore_grnd_img-fore_grnd_img_inner_pointwise_lowest;

% erosion of small blobs by opening
fore_grnd_img_inner_markers_open = fore_grnd_img_inner_markers_pointwise_lowest; %imopen(fore_grnd_img_inner_markers_pointwise_lowest,strel_small_ring);

% inner markers with original values
fore_grnd_img_inner_markers = fore_grnd_img.* (fore_grnd_img_inner_markers_open > 0);

% blob labeling of inner markers
[blob_label_inner_marker_image,blob_count_inner_marker] = bwlabel(fore_grnd_img_inner_markers,4);

% %%%%%
% figure(100);
% imagesc(blob_label_inner_marker_image), colormap(gray);
% set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
% set(gca, 'PlotBoxAspectRatio',[1 1 1])
% set(gcf,'Position',[-1854 645 550 418]);
% figure;
% imagesc(fore_grnd_img_inner_markers), colormap(gray);
% set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
% set(gca, 'PlotBoxAspectRatio',[1 1 1])
% set(gcf,'Position',[-1480 678 550 418]);
% %%%%%


% Startup:
% find max intensity and its coordinates
while any(blob_label_inner_marker_image(:))
    % find max intensity and its coordinates
    [max_int, max_index] = max(fore_grnd_img_inner_markers(:));
    [max_y, max_x] = ind2sub(size(fore_grnd_img), max_index);
    
    % Update window coordinates
    MB_window_coord(1,1) = max_y-MB_window_size_search_new(1);
    MB_window_coord(2,1) = max_y+MB_window_size_search_new(1);
    MB_window_coord(1,2) = max_x-MB_window_size_search_new(2);
    MB_window_coord(2,2) = max_x+MB_window_size_search_new(2);
    
    
    if ( ~((MB_window_coord(1,1) <= 0) && (MB_window_coord(2,1) > size(fore_grnd_img,1)) && (MB_window_coord(1,2) <= 0) && (MB_window_coord(2,2) > size(fore_grnd_img,2)) ))
    % check for out of bounce
    % y1
    if MB_window_coord(1,1) <= 0
        MB_window_coord(1,1) = 1;
    end
    % y2
    if MB_window_coord(2,1) > size(fore_grnd_img,1)
        MB_window_coord(2,1) = size(fore_grnd_img,1);
    end
    % x1
    if MB_window_coord(1,2) <= 0
        MB_window_coord(1,2) = 1;
    end
    % x2
    if MB_window_coord(2,2) > size(fore_grnd_img,2)
        MB_window_coord(2,2) = size(fore_grnd_img,2);
    end
    
    % Create temp window
    fore_grnd_img_temp_window = fore_grnd_img(MB_window_coord(1,1):MB_window_coord(2,1),MB_window_coord(1,2):MB_window_coord(2,2));
    
    % FWHM thresholding
    fore_grnd_img_temp_window(find(fore_grnd_img_temp_window < max_int/3)) = 0;

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
    
    % Count number of blobs in FWHM_window
    [blob_label_foregrnd_image_window,blob_count_foregrnd_window] = bwlabel(fore_grnd_img_temp_window,4);
    % Find label value for max intensity point
    blob_label_value = blob_label_foregrnd_image_window(find(fore_grnd_img_temp_window == max(fore_grnd_img_temp_window(:))));
    % remove all other blobs than max intensity
    blob_label_foregrnd_image_window(find(blob_label_foregrnd_image_window ~= blob_label_value(1)))  = 0; 
        
    % Check for MB_count condition
    if (blob_count_inner_marker_window <= MB_count_condition) && (blob_count_foregrnd_window <= MB_count_condition)
        
        % Features BW: Area, Centroid, BoundingBox, Eccentricity, Extent,
        % FilledArea, MajorAxisLength, MinorAxisLength, Orientatiom, Perimeter 
        MB_blob_features_BW = regionprops(blob_label_foregrnd_image_window,...
            'Area','BoundingBox','Centroid','Eccentricity','Extent','MajorAxisLength','MinorAxisLength','Orientation','Perimeter');

        % Features I: MaxIntensity, MeanIntensity, WeigtedCentroid
        MB_blob_features_I = regionprops(blob_label_foregrnd_image_window,fore_grnd_img_temp_window,'MaxIntensity','MeanIntensity','WeightedCentroid');

        % updata MB
        MB(MB_index).state = 1;
        MB(MB_index).old_pos = [0, 0];
        MB(MB_index).new_pos = [max_y, max_x];
        MB(MB_index).vel = [0,0];
        MB(MB_index).id = MB_index;
        MB(MB_index).max_int = MB_blob_features_I.MaxIntensity;
        MB(MB_index).mean_int = MB_blob_features_I.MeanIntensity;
        MB(MB_index).centroid_I = fliplr(round(MB_blob_features_I.WeightedCentroid)) + [MB_window_coord(1,1)-1, MB_window_coord(1,2)-1];
        MB(MB_index).centroid_BW = fliplr(round(MB_blob_features_BW.Centroid)) + [MB_window_coord(1,1)-1, MB_window_coord(1,2)-1];
        MB(MB_index).area = MB_blob_features_BW.Area;
        MB(MB_index).bounding_box = MB_blob_features_BW.BoundingBox;
        MB(MB_index).eccentricity = MB_blob_features_BW.Eccentricity;
        MB(MB_index).extent = MB_blob_features_BW.Extent;
        MB(MB_index).major_ax = MB_blob_features_BW.MajorAxisLength;
        MB(MB_index).minor_ax = MB_blob_features_BW.MinorAxisLength;
        MB(MB_index).orientation = MB_blob_features_BW.Orientation;
        MB(MB_index).perimeter = MB_blob_features_BW.Perimeter;
        MB(MB_index).count = [blob_count_inner_marker_window blob_count_foregrnd_window];
        MB(MB_index).age = [idx_seq idx_seq 1];

        MB_index = MB_index + 1;
        %figure, imagesc(fore_grnd_img_temp(:,:)), colormap(gray);
    end
end

% make a copy of MB for later logging
MB_log = MB;
% after startup

steps = 1; % number of frames skipped
% Number of frame steps set
for idx_seq = idx_seq_start + steps:steps:150%size(fore_grnd_img,3)
    % New frame loaded
    fore_grnd_img = load_img_tracking_MB_res1010(idx_seq);
    
    % inner maker image
    % dilate
    fore_grnd_img_inner_marker_dilate = imdilate(fore_grnd_img,strel_ring);
 
    % logical pointwise lowest identifier
    pointwise_mask = (fore_grnd_img_inner_marker_dilate > fore_grnd_img);
   
    % pointwise lowest
    fore_grnd_img_inner_pointwise_lowest = fore_grnd_img.*pointwise_mask + fore_grnd_img_inner_marker_dilate.*(~pointwise_mask); 
   
    % (original - pointwise lowest)
    fore_grnd_img_inner_markers_pointwise_lowest = fore_grnd_img-fore_grnd_img_inner_pointwise_lowest;
    
    % erosion of small blobs by opening
    fore_grnd_img_inner_markers_open = fore_grnd_img_inner_markers_pointwise_lowest; % imopen(fore_grnd_img_inner_markers_pointwise_lowest,strel_small_ring);
    
    % inner markers with original values
    fore_grnd_img_inner_markers = fore_grnd_img.* (fore_grnd_img_inner_markers_open > 0);
    
    % blob labeling of inner markers
    [blob_label_inner_marker_image,blob_count_inner_marker] = bwlabel(fore_grnd_img_inner_markers,4);
        
    % sort MB's after intensity
    MB_sort_index = logical([MB.state]);
    [~,MB_sort_index] = sort([MB.max_int].*MB_sort_index,'descend'); 
    
    % Search for movement of exiting MB's
    for i = 1:sum([MB.state]) 
        % update MB_index
        MB_index = MB_sort_index(i);
        
        % Update window coordinates
        MB_window_coord(1,1) = MB(MB_index).new_pos(1)-MB_window_size_search_existing(1)+MB(MB_index).vel(1);
        MB_window_coord(2,1) = MB(MB_index).new_pos(1)+MB_window_size_search_existing(1)+MB(MB_index).vel(1);
        MB_window_coord(1,2) = MB(MB_index).new_pos(2)-MB_window_size_search_existing(2)+MB(MB_index).vel(2);
        MB_window_coord(2,2) = MB(MB_index).new_pos(2)+MB_window_size_search_existing(2)+MB(MB_index).vel(2);

        % check for out of bounce
        % y1
        if MB_window_coord(1,1) <= 0
            MB_window_coord(1,1) = 1;
        end
        % y2
        if MB_window_coord(2,1) > size(fore_grnd_img,1)
            MB_window_coord(2,1) = size(fore_grnd_img,1);
        end
        % x1
        if MB_window_coord(1,2) <= 0
            MB_window_coord(1,2) = 1;
        end
        % x2
        if MB_window_coord(2,2) > size(fore_grnd_img,2)
            MB_window_coord(2,2) = size(fore_grnd_img,2);
        end
        
        % Get blob window
        blob_label_image_inner_marker_window = blob_label_inner_marker_image(MB_window_coord(1,1):MB_window_coord(2,1),MB_window_coord(1,2):MB_window_coord(2,2));        

        if any(blob_label_image_inner_marker_window(:))
            
            % Count number of blobs in inner_marker_window
            blobs_in_window = unique(blob_label_image_inner_marker_window);
            blobs_in_window(1) = [];
            blob_count_inner_marker_window = length(blobs_in_window);
        
            % Create temp window around the respective MB
            fore_grnd_img_temp_window = fore_grnd_img(MB_window_coord(1,1):MB_window_coord(2,1),MB_window_coord(1,2):MB_window_coord(2,2));

            % Create temp window around the respective MB from inner markers
            fore_grnd_img_inner_markers_window = fore_grnd_img_inner_markers(MB_window_coord(1,1):MB_window_coord(2,1),MB_window_coord(1,2):MB_window_coord(2,2));
            
            [max_int, max_index_window] = max(fore_grnd_img_inner_markers_window(:));
            [max_y_window, max_x_window] = ind2sub(size(fore_grnd_img_temp_window), max_index_window);
            max_y = max_y_window + MB_window_coord(1,1)-1;
            max_x = max_x_window + MB_window_coord(1,2)-1;
            
            % FWHM thresholding
            fore_grnd_img_temp_window(find(fore_grnd_img_temp_window < max_int/threshold)) = 0;
            
            % Picks out the blobs inside the blob_label_inner_maker_window
            for i = 1:blob_count_inner_marker_window
                fore_grnd_img_inner_markers(find(blob_label_inner_marker_image == blobs_in_window(i))) = 0;
                blob_label_inner_marker_image(find(blob_label_inner_marker_image == blobs_in_window(i))) = 0;
            end
            
            % Count number of blobs in FWHM_window
            [blob_label_foregrnd_image_window,blob_count_foregrnd_window] = bwlabel(fore_grnd_img_temp_window,4);
            % Find label value for max intensity point
            blob_label_value = blob_label_foregrnd_image_window(find(fore_grnd_img_temp_window == max(fore_grnd_img_temp_window(:))));
            % remove all other blobs than max intensity
            blob_label_foregrnd_image_window(find(blob_label_foregrnd_image_window ~= blob_label_value(1))) = 0; 


            if (blob_count_inner_marker_window <= MB_count_condition) && (blob_count_foregrnd_window <= MB_count_condition)

                % Features BW: Area, Centroid, BoundingBox, Eccentricity, Extent,
                % FilledArea, MajorAxisLength, MinorAxisLength, Orientatiom, Perimeter 
                MB_blob_features_BW = regionprops(blob_label_foregrnd_image_window,...
                    'Area','BoundingBox','Centroid','Eccentricity','Extent','MajorAxisLength','MinorAxisLength','Orientation','Perimeter');

                % Features I: MaxIntensity, MeanIntensity, WeigtedCentroid
                MB_blob_features_I = regionprops(blob_label_foregrnd_image_window, fore_grnd_img_temp_window,'MaxIntensity','MeanIntensity','WeightedCentroid');

                % updata MB
                MB(MB_index).old_pos = MB(MB_index).new_pos;
                MB(MB_index).new_pos = [max_y,max_x];
                MB(MB_index).vel = MB(MB_index).new_pos-MB(MB_index).old_pos;
                MB(MB_index).max_int = MB_blob_features_I.MaxIntensity;
                MB(MB_index).mean_int = MB_blob_features_I.MeanIntensity;
                MB(MB_index).centroid_I = fliplr(round(MB_blob_features_I.WeightedCentroid)) + [MB_window_coord(1,1)-1, MB_window_coord(1,2)-1];
                MB(MB_index).centroid_BW = fliplr(round(MB_blob_features_BW.Centroid)) + [MB_window_coord(1,1)-1, MB_window_coord(1,2)-1];
                MB(MB_index).area = MB_blob_features_BW.Area;
                MB(MB_index).bounding_box = MB_blob_features_BW.BoundingBox;
                MB(MB_index).eccentricity = MB_blob_features_BW.Eccentricity;
                MB(MB_index).extent = MB_blob_features_BW.Extent;
                MB(MB_index).major_ax = MB_blob_features_BW.MajorAxisLength;
                MB(MB_index).minor_ax = MB_blob_features_BW.MinorAxisLength;
                MB(MB_index).orientation = MB_blob_features_BW.Orientation;
                MB(MB_index).perimeter = MB_blob_features_BW.Perimeter;
                MB(MB_index).count = [blob_count_inner_marker_window blob_count_foregrnd_window];
                MB(MB_index).age = MB(MB_index).age + [0 steps 1];

                % updata MB_log
                MB_log(MB_index).old_pos(MB(MB_index).age(3),:) = MB(MB_index).old_pos;
                MB_log(MB_index).new_pos(MB(MB_index).age(3),:) = MB(MB_index).new_pos;
                MB_log(MB_index).vel(MB(MB_index).age(3),:) = MB(MB_index).vel;
                MB_log(MB_index).max_int(MB(MB_index).age(3)) = MB(MB_index).max_int;
                MB_log(MB_index).mean_int(MB(MB_index).age(3)) = MB(MB_index).mean_int;
                MB_log(MB_index).centroid_I(MB(MB_index).age(3),:) = MB(MB_index).centroid_I;
                MB_log(MB_index).centroid_BW(MB(MB_index).age(3),:) = MB(MB_index).centroid_BW;
                MB_log(MB_index).area(MB(MB_index).age(3)) = MB(MB_index).area;
                MB_log(MB_index).bounding_box(MB(MB_index).age(3),:) = MB(MB_index).bounding_box;
                MB_log(MB_index).eccentricity(MB(MB_index).age(3)) = MB(MB_index).eccentricity;
                MB_log(MB_index).extent(MB(MB_index).age(3)) = MB(MB_index).extent;
                MB_log(MB_index).major_ax(MB(MB_index).age(3)) = MB(MB_index).major_ax;
                MB_log(MB_index).minor_ax(MB(MB_index).age(3)) = MB(MB_index).minor_ax;
                MB_log(MB_index).orientation(MB(MB_index).age(3)) = MB(MB_index).orientation;
                MB_log(MB_index).perimeter(MB(MB_index).age(3)) = MB(MB_index).perimeter;
                MB_log(MB_index).count(MB(MB_index).age(3),:) = MB(MB_index).count;
                MB_log(MB_index).age = MB(MB_index).age;
            else
                % Remove MB
                MB_log(MB_index).state = 0;
            end
        else
            % Remove MB
            MB_log(MB_index).state = 0;
        end
        % Check next MB
    end
    
    % Seach for new MB's
    while any(blob_label_inner_marker_image(:))         
        % find max intensity and its coordinates
        [max_int, max_index] = max(fore_grnd_img_inner_markers(:));
        [max_y, max_x] = ind2sub(size(fore_grnd_img), max_index);

        % Update window coordinates
        MB_window_coord(1,1) = max_y-MB_window_size_search_new(1);
        MB_window_coord(2,1) = max_y+MB_window_size_search_new(1);
        MB_window_coord(1,2) = max_x-MB_window_size_search_new(2);
        MB_window_coord(2,2) = max_x+MB_window_size_search_new(2);

        % check for out of bounce
        % y1
        if MB_window_coord(1,1) <= 0
            MB_window_coord(1,1) = 1;
        end
        % y2
        if MB_window_coord(2,1) > size(fore_grnd_img,1)
            MB_window_coord(2,1) = size(fore_grnd_img,1);
        end
        % x1
        if MB_window_coord(1,2) <= 0
            MB_window_coord(1,2) = 1;
        end
        % x2
        if MB_window_coord(2,2) > size(fore_grnd_img,2)
            MB_window_coord(2,2) = size(fore_grnd_img,2);
        end

        % Create temp window
        fore_grnd_img_temp_window = fore_grnd_img(MB_window_coord(1,1):MB_window_coord(2,1),MB_window_coord(1,2):MB_window_coord(2,2));

        % FWHM thresholding
        fore_grnd_img_temp_window(find(fore_grnd_img_temp_window < max_int/threshold)) = 0;

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

        % Count number of blobs in FWHM_window
        [blob_label_foregrnd_image_window,blob_count_foregrnd_window] = bwlabel(fore_grnd_img_temp_window,4);
        % Find label value for max intensity point
        blob_label_value = blob_label_foregrnd_image_window(find(fore_grnd_img_temp_window == max(fore_grnd_img_temp_window(:))));
        % remove all other blobs than max intensity
        blob_label_foregrnd_image_window(find(blob_label_foregrnd_image_window ~= blob_label_value(1)))  = 0; 
        
        % Check for MB_count condition
        if (blob_count_inner_marker_window <= MB_count_condition) && (blob_count_foregrnd_window <= MB_count_condition)
            % Index update
            MB_index = max([MB_log(:).id]) + 1;
            
            % Picks out the blob with max intensity
            blob_label_image_extract_blobs = zeros(size(blob_label_inner_marker_image));
            blob_label_image_extract_blobs(find(blob_label_inner_marker_image == blob_label_inner_marker_image(max_y,max_x))) = 1; % max_y, max_x legit here..

            % Features BW: Area, Centroid, BoundingBox, Eccentricity, Extent,
            % FilledArea, MajorAxisLength, MinorAxisLength, Orientatiom, Perimeter 
            MB_blob_features_BW = regionprops(blob_label_foregrnd_image_window,...
                'Area','BoundingBox','Centroid','Eccentricity','Extent','MajorAxisLength','MinorAxisLength','Orientation','Perimeter');

            % Features I: MaxIntensity, MeanIntensity, WeigtedCentroid
            MB_blob_features_I = regionprops(blob_label_foregrnd_image_window,fore_grnd_img_temp_window,'MaxIntensity','MeanIntensity','WeightedCentroid');

            % Set respective BLOB to zero
            fore_grnd_img = fore_grnd_img.*(~blob_label_image_extract_blobs);

            % updata MB
            MB(MB_index).state = 1;
            MB(MB_index).old_pos = [0, 0];
            MB(MB_index).new_pos = [max_y, max_x];
            MB(MB_index).vel = [0,0];
            MB(MB_index).id = MB_index;
            MB(MB_index).max_int = MB_blob_features_I.MaxIntensity;
            MB(MB_index).mean_int = MB_blob_features_I.MeanIntensity;
            MB(MB_index).centroid_I = fliplr(round(MB_blob_features_I.WeightedCentroid)) + [MB_window_coord(1,1)-1, MB_window_coord(1,2)-1];
            MB(MB_index).centroid_BW = fliplr(round(MB_blob_features_BW.Centroid)) + [MB_window_coord(1,1)-1, MB_window_coord(1,2)-1];
            MB(MB_index).area = MB_blob_features_BW.Area;
            MB(MB_index).bounding_box = MB_blob_features_BW.BoundingBox;
            MB(MB_index).eccentricity = MB_blob_features_BW.Eccentricity;
            MB(MB_index).extent = MB_blob_features_BW.Extent;
            MB(MB_index).major_ax = MB_blob_features_BW.MajorAxisLength;
            MB(MB_index).minor_ax = MB_blob_features_BW.MinorAxisLength;
            MB(MB_index).orientation = MB_blob_features_BW.Orientation;
            MB(MB_index).perimeter = MB_blob_features_BW.Perimeter;
            MB(MB_index).count = [blob_count_inner_marker_window blob_count_foregrnd_window];
            MB(MB_index).age = [idx_seq idx_seq 1];

            % updata MB_log
            MB_log(MB_index).state = 1;
            MB_log(MB_index).old_pos = MB(MB_index).old_pos;
            MB_log(MB_index).new_pos = MB(MB_index).new_pos;
            MB_log(MB_index).vel = MB(MB_index).vel;
            MB_log(MB_index).id = MB_index;
            MB_log(MB_index).max_int = MB(MB_index).max_int;
            MB_log(MB_index).mean_int = MB(MB_index).mean_int;
            MB_log(MB_index).centroid_I = MB(MB_index).centroid_I;
            MB_log(MB_index).centroid_BW = MB(MB_index).centroid_BW;
            MB_log(MB_index).area = MB(MB_index).area;
            MB_log(MB_index).bounding_box = MB(MB_index).bounding_box;
            MB_log(MB_index).eccentricity = MB(MB_index).eccentricity;
            MB_log(MB_index).extent = MB(MB_index).extent;
            MB_log(MB_index).major_ax = MB(MB_index).major_ax;
            MB_log(MB_index).minor_ax = MB(MB_index).minor_ax;
            MB_log(MB_index).orientation = MB(MB_index).orientation;
            MB_log(MB_index).perimeter = MB(MB_index).perimeter;
            MB_log(MB_index).count = MB(MB_index).count;
            MB_log(MB_index).age = MB(MB_index).age;
        end
    end
end

% Running time
total_time = toc % total time

% remove first instances of pos and vel due to startup
for MB_index = 1:length(MB)
    MB_log(MB_index).old_pos(1,:) = [];
    MB_log(MB_index).vel(1,:) = [];
end

% Load_img parameters:
n_bck_grnd = 5; % # of img's to make back ground
n_bck_grnd_skip = 10; % # of img's skipped between each background img

% save MB_log with variable used
% get path to current folder
[current_dir,current_file,current_ext] = fileparts([mfilename('fullpath') '.m']) ;
 %%
% check for existing file
for i = 0:10
    file_name = [current_dir '/' current_file '_' num2str(i,'%d') '.mat'];
    file_exist_check = exist(file_name);

    if file_exist_check == 0
        fprintf(['Data saved in: ' file_name '\n']);
        save(file_name, 'MB_log','MB_count_condition','MB_age_condition','MB_dens_condition_avg','MB_dens_condition_single','MB_window_size_search_new','MB_window_size_search_existing','img_size','steps','n_bck_grnd','n_bck_grnd_skip')
        break;
    end
end

%              figure;
%             imagesc(blob_label_inner_marker_image), colormap(gray);
%             set(gca, 'DataAspectRatio',[1 1 1]) % set data aspect ratio in zoom box
%             set(gca, 'PlotBoxAspectRatio',[1 1 1])
%             set(gcf,'Position',[-1854 645 550 418]);
