log_compression = @(x) 20*log10(x/max(x(:)));

% Folder path to data 
folderName = '/data/cfudata3/ramosh/cavh/microbubble-experiments/Part6-PM_Half_Front_10V_6_5MHz_2Beam1Line_Flow_Ramp_5uls_to_0uls_180s_1to100';
% folderName = '/data/cfudata3/ramosh/cavh/microbubble-experiments/twochannel_5ul_sec';



%% Read usecase
usecase_xml = xmlread([folderName '/usecase.xml']);
% read sampling rate
temp = usecase_xml.getElementsByTagName('receiveSampleFrequency');
sampling_rate = str2double(temp.item(0).getFirstChild.getData);
% read start line
temp = usecase_xml.getElementsByTagName('actualStartLine');
start_line = str2double(temp.item(0).getFirstChild.getData);
% read end line
temp = usecase_xml.getElementsByTagName('actualStopLine');
end_line = str2double(temp.item(0).getFirstChild.getData);
% read element pitch
temp = usecase_xml.getElementsByTagName('pitch');
element_pitch = str2double(temp.item(0).getFirstChild.getData);
% read line density
temp = usecase_xml.getElementsByTagName('lineDensity');
line_density = str2double(temp.item(0).getFirstChild.getData);
% read stop depth
temp = usecase_xml.getElementsByTagName('actualStopDepth');
stop_depth = str2double(temp.item(0).getFirstChild.getData);
% read wave sample to calculate center_frequency
Tx_transmit_rate = 120e6;
pulse_per_wawe = 1.5;
temp = usecase_xml.getElementsByTagName('sampleCount');
center_frequency = (str2double(temp.item(0).getFirstChild.getData)/(Tx_transmit_rate*pulse_per_wawe))^(-1);

%% Image parameters

% Area of interst
area_of_interest.axial_init= 950; 
area_of_interest.axial_end= 1250; 
area_of_interest.lateral_init= 1; 
area_of_interest.lateral_end= 60; 

% area_of_interest.axial_init= 775; 
% area_of_interest.axial_end= 950; 
% area_of_interest.lateral_init= 1; 
% area_of_interest.lateral_end= 60; 

% area_of_interest.axial_init= 700; 
% area_of_interest.axial_end= 1000; 
% area_of_interest.lateral_init= 1; 
% area_of_interest.lateral_end= 60; 

% Wanted resolution
img_resolution.axial_new = 10e-6;
img_resolution.lateral_new = 10e-6;

img_type = 'PI';


%% Calculate frame rate
[~, t1_temp, ~] = load_img(1, img_type, folderName, img_resolution, area_of_interest);
[~, t2_temp, ~] = load_img(2, img_type, folderName, img_resolution, area_of_interest);
fps = 1/(t2_temp/(1e6)-t1_temp/(1e6));

%% Get img size

img_size = size(load_img(1, img_type, folderName, img_resolution, area_of_interest));

%% Tracking parameters
% Tracking parameters
MB_window_coord_search = zeros(2,2); % Coordinates for search window
MB_window_coord_search_clear = zeros(2,2); % Coordinates for search clear window

MB_window_size_search_localization = [25 25]; % [y,x] Search window for localization of PSF's
MB_window_size_search_new = [20 20]; % [y,x] Search window for new MB's
MB_window_size_search_existing = [15 15]; % [y,x] Search window for ''old'' MB's
MB_window_size_search_clear = [30 30]; % [y,x] Search window for ''old'' MB's
MB_window_threshold = 1.3; % Window Threshold (actual MB_window_threshold = Max_intensity*1/threshold)


idx_frame_start = 500;
n_frame = 15000;

% initializing MB struct

MB_log = []; % Current MB
MB_log.state = 0;
MB_log.old_pos = zeros(1,2);
MB_log.new_pos = zeros(1,2);
MB_log.pred_pos = zeros(1,2);
MB_log.vel = zeros(1,2); 
MB_log.max_int = 0;
MB_log.centroid = zeros(1,2);
MB_log.area = 0;
MB_log.eccentricity = 0;
MB_log.orientation = 0;
MB_log.perimeter = 0;
MB_log.id = 0;
MB_log.count = 0;
% MB_log = MB;

MB_idx = 0; % Index for each MB
MB_window_out_of_bounce_flag = 0; % Checks if search windows is outside image

tic
for idx_frame=idx_frame_start:idx_frame_start+n_frame
    idx_frame
    
    % Load img
    [PI_img, time_stamp, line_count] = load_img(idx_frame, img_type, folderName, img_resolution, area_of_interest);

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
    [~,MB_sort_index] = sort([MB_log.state],'descend');
  
    % Connect MB's
    for i = 1:sum(logical([MB_log.state]))
        % Update MB_index
        MB_idx = MB_sort_index(i);
        
        % Last recorded position:
        MB_pos_y = MB_log(MB_idx).new_pos(end,1);
        MB_pos_x = MB_log(MB_idx).new_pos(end,2);
        
        % Update window coordinates
        if MB_log(MB_idx).age(3) == 1
            MB_window_coord_search(1,1) = MB_pos_y-MB_window_size_search_new(1);
            MB_window_coord_search(2,1) = MB_pos_y+MB_window_size_search_new(1);
            MB_window_coord_search(1,2) = MB_pos_x-MB_window_size_search_new(2);
            MB_window_coord_search(2,2) = MB_pos_x+MB_window_size_search_new(2);
        else
            MB_window_coord_search(1,1) = MB_pos_y-MB_window_size_search_existing(1)+MB_log(MB_idx).vel(1);
            MB_window_coord_search(2,1) = MB_pos_y+MB_window_size_search_existing(1)+MB_log(MB_idx).vel(1);
            MB_window_coord_search(1,2) = MB_pos_x-MB_window_size_search_existing(2)+MB_log(MB_idx).vel(2);
            MB_window_coord_search(2,2) = MB_pos_x+MB_window_size_search_existing(2)+MB_log(MB_idx).vel(2);
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
%       figure(); imagesc(img_temp_window); colormap('gray');%---
        
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
            
%------------------------------------------------
            % Load img
%             temp_PI_img = load_img(idx_frame-1, img_type, folderName, img_resolution, area_of_interest);
%             temp_img_temp_window = temp_PI_img(MB_window_coord_search(1,1):MB_window_coord_search(2,1),MB_window_coord_search(1,2):MB_window_coord_search(2,2));
%             figure(); imagesc(temp_img_temp_window); colormap('gray');%---
%-------------------------------------------------


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
                
                % Updata MB log
                MB_log(MB_idx).age = MB_log(MB_idx).age + [0 1 1];
                MB_log(MB_idx).state = max_int;
                MB_log(MB_idx).old_pos(MB_log(MB_idx).age(3),:) = MB_log(MB_idx).new_pos(end,:);
                MB_log(MB_idx).new_pos(MB_log(MB_idx).age(3),:) = fliplr(round(MB_blob_features.WeightedCentroid)) + [MB_window_coord_search(1,1)-1, MB_window_coord_search(1,2)-1]; %[max_y,max_x];
                MB_log(MB_idx).pred_pos(MB_log(MB_idx).age(3),:) = [MB_pos_y+MB_log(MB_idx).vel(1),MB_pos_x+MB_log(MB_idx).vel(2);];
                MB_log(MB_idx).vel(MB_log(MB_idx).age(3),:) = MB_log(MB_idx).new_pos(end,:)-MB_log(MB_idx).old_pos(end,:);
                MB_log(MB_idx).max_int(MB_log(MB_idx).age(3)) = max_int;
                MB_log(MB_idx).centroid(MB_log(MB_idx).age(3),:) = fliplr(round(MB_blob_features.WeightedCentroid)) + [MB_window_coord_search(1,1)-1, MB_window_coord_search(1,2)-1];
                MB_log(MB_idx).area(MB_log(MB_idx).age(3)) = MB_blob_features_bw.Area;
                MB_log(MB_idx).eccentricity(MB_log(MB_idx).age(3)) = MB_blob_features_bw.Eccentricity;
                MB_log(MB_idx).orientation(MB_log(MB_idx).age(3)) = MB_blob_features_bw.Orientation;
                MB_log(MB_idx).perimeter(MB_log(MB_idx).age(3)) = MB_blob_features_bw.Perimeter;
                MB_log(MB_idx).count(MB_log(MB_idx).age(3),:) = [blob_count_global_window blob_count_window];
                
            else
                % Remove MB from list
                MB_log(MB_idx).state = 0;
            end
        else
            % Remove MB from list
            MB_log(MB_idx).state = 0;
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
            MB_idx = max([MB_log(:).id]) + 1;
            
            % Features: MaxIntensity, MeanIntensity, WeigtedCentroid
            MB_blob_features = regionprops(img_blob_label_window, img_temp_window,'MaxIntensity','WeightedCentroid');
            
            % Features: Area, Eccentricity, Orientation, Perimeter
            MB_blob_features_bw = regionprops(img_blob_label_window,'Area','Eccentricity','Orientation','Perimeter');
            
            % Updata MB log
            MB_log(MB_idx).state = max_int;
            MB_log(MB_idx).old_pos = [0, 0];
            MB_log(MB_idx).new_pos = fliplr(round(MB_blob_features.WeightedCentroid)) + [MB_window_coord_search(1,1)-1, MB_window_coord_search(1,2)-1];%[max_y, max_x];
            MB_log(MB_idx).pred_pos = [0,0];
            MB_log(MB_idx).vel = [0,0];
            MB_log(MB_idx).id = MB_idx;
            MB_log(MB_idx).max_int = max_int;
            MB_log(MB_idx).centroid = fliplr(round(MB_blob_features.WeightedCentroid)) + [MB_window_coord_search(1,1)-1, MB_window_coord_search(1,2)-1];
            MB_log(MB_idx).area = MB_blob_features_bw.Area;
            MB_log(MB_idx).eccentricity = MB_blob_features_bw.Eccentricity;
            MB_log(MB_idx).orientation = MB_blob_features_bw.Orientation;
            MB_log(MB_idx).perimeter = MB_blob_features_bw.Perimeter;
            MB_log(MB_idx).count = [blob_count_global_window blob_count_window];
            MB_log(MB_idx).age = [idx_frame idx_frame 1];
        end
    end
end

% Remove first instances of pos and vel due to startup
for MB_idx = 1:length(MB_log)
    MB_log(MB_idx).old_pos(1,:) = [];
    MB_log(MB_idx).vel(1,:) = [];
end

% Save all parameters
MB_data.folderName = folderName;
MB_data.MB_log = MB_log;
MB_data.samplingRate = sampling_rate;
MB_data.startLine = start_line;
MB_data.endLine = end_line;
MB_data.elementPitch = element_pitch;
MB_data.lineDensity = line_density;
MB_data.stopDepth = stop_depth;
MB_data.txTransmitFrequency = Tx_transmit_rate;
MB_data.centerFrequency = center_frequency;
MB_data.fps = fps;
MB_data.imgResolution = img_resolution;
MB_data.imgSize = img_size;
MB_data.imgType = img_type;
MB_data.areaOfInterest = area_of_interest;
MB_data.MB_window_size_search_localization = MB_window_size_search_localization;
MB_data.MB_window_size_search_new = MB_window_size_search_new;
MB_data.MB_window_size_search_existing = MB_window_size_search_existing;
MB_data.MB_window_size_search_clear = MB_window_size_search_clear;
MB_data.MB_window_threshold = MB_window_threshold;
MB_data.startFrame = idx_frame_start;
MB_data.numberOfFrames = n_frame;


toc

% Save parameters in MB_tracking_parameter

%MB_tracking_parameter.MB_log = MB_log;



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


figure();
% subplot(1,2,1)
PI_img = load_img(5000, img_type, folderName, img_resolution, area_of_interest);
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
figure(7); 
PI_img = load_img(5000, img_type, folderName, img_resolution, area_of_interest);;


outputVideo=VideoWriter('test5_vid');
outputVideo.FrameRate=3;
open(outputVideo);
mov(1:50)= struct('cdata',[],'colormap',[]);
pause(0.3)
frame_timestamp_linecount = [];
idx = 0;
img_type = 'PI'
for idx_frame = 15000:22000
    idx_frame
    idx = idx + 1;
    [PI_img, time_stamp, line_count] = load_img(idx_frame, img_type, folderName, img_resolution, area_of_interest);
%     imagesc(log_compression(abs(PI_img)),[-60 0])
    frame_timestamp_linecount(idx,:) = [time_stamp, line_count, idx_frame];
    if idx>1
        1/(diff(frame_timestamp_linecount(idx-1:idx,1),1,1)*1e-6);
    end
    
    figure(7);
    imagesc(lateral_axis*1000,depth_axis*1000,log_compression(abs(PI_img)),[-20 0])
%     imagesc(lateral_axis*1000,depth_axis*1000,abs(PI_img),[0 2000])
    
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






