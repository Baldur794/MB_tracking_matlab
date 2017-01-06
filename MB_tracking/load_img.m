function [ img, new_time_stamp, new_line_count ] = load_img(idx_frame, img_type, folderName, img_resolution, area_of_interest)
% Outputs foreground img.

addpath('/home/cavh/Documents/papers/Journals/SuperResolution/scripts/Reading_bk3000/library/BkmDataFile');

% Depends on number of images per frame
streamName{1}='export-0.bkmdf';
streamName{2}='export-1.bkmdf';
% streamName{3}='export-2.bkmdf';
% streamName{4}='export-3.bkmdf';

% Read usecase
usecase_xml = xmlread([folderName '/usecase.xml']);
% read sampling rate
temp = usecase_xml.getElementsByTagName('receiveSampleFrequency');
sampling_rate = str2double(temp.item(0).getFirstChild.getData);
% read element pitch
temp = usecase_xml.getElementsByTagName('pitch');
element_pitch = str2double(temp.item(0).getFirstChild.getData);
% read line density
temp = usecase_xml.getElementsByTagName('lineDensity');
line_density = str2double(temp.item(0).getFirstChild.getData);

% Propagation speed of sound
c = 1540;

% Original resolution
img_resolution.axial_original = 1/sampling_rate*c/2;
img_resolution.lateral_original = element_pitch/line_density;

% To obtain wanted resolution
interpolate_flag = true;
if interpolate_flag
    interpolation_factor.axial = img_resolution.axial_original/img_resolution.axial_new;
    interpolation_factor.lateral = img_resolution.lateral_original/img_resolution.lateral_new;
else
    interpolation_factor_axial = 1;
    interpolation_factor_lateral = 1;
end
interpolation_type = 'spline'; % Interpolation method

% Grid for interpolation
[X,Y] = meshgrid(area_of_interest.lateral_init:area_of_interest.lateral_end,area_of_interest.axial_init:area_of_interest.axial_end);
[Xq,Yq] = meshgrid(area_of_interest.lateral_init:1/interpolation_factor.lateral:area_of_interest.lateral_end,area_of_interest.axial_init:1/interpolation_factor.axial:area_of_interest.axial_end);

% Filter parameters
n_bck_grnd_skip = 100; % Distance between frame in interest and bck image
n_bck_jump = 20; % Number of frames between each img used for final bck image
n_bck = 5; % Number of imgs used for bck image
n_fore = 1; % Number of imgs used for fore image

% Bandpass filter in RF domain
bpf.filt_order=100;            % Order of filter
bpf.cuttoff1= 2e6;             % Lower cutoff MHz
bpf.cuttoff2= 18e6;            % Upper cutoff MHz
bpf.filter = designfilt('bandpassfir', 'FilterOrder', bpf.filt_order, 'CutoffFrequency1', bpf.cuttoff1, 'CutoffFrequency2', bpf.cuttoff2, 'SampleRate', sampling_rate);

if strcmp('PI',img_type) % For PI image
    % Choose which images to use for PI: [x y] -> PI = x-y
    select_img = [1 2];
    
    % % Background img
    bck_grnd_img = zeros(size(Xq));
    
    if n_bck > 0
        % sum of n_bck_grnd images
        for it_frame = idx_frame-n_bck_grnd_skip-(n_bck-1)*n_bck_jump:n_bck_jump:idx_frame-n_bck_grnd_skip
            % Ensure valid index
            if it_frame <= 0
                it_frame = 1;
            end
            for it_stream=1:2
                [header temp]=bkmload([folderName filesep streamName{select_img(it_stream)}],it_frame);
                temp_data{it_stream} = hilbert(double(temp(area_of_interest.axial_init:area_of_interest.axial_end,area_of_interest.lateral_init:area_of_interest.lateral_end)));
                temp_data{it_stream} = filtfilt(bpf.filter,temp_data{it_stream});
            end
            PI_img_bck = abs(temp_data{select_img(1)}+temp_data{select_img(2)});
            PI_img_bck = interp2(X,Y,PI_img_bck,Xq,Yq,interpolation_type);
            
            bck_grnd_img = bck_grnd_img + PI_img_bck;
        end
        % background average
        bck_grnd_img = bck_grnd_img/n_bck;
    end
    
    % Foreground img
    fore_grnd_img = zeros(size(Xq));
    
    % Sum of n_fore images
    for it_frame = idx_frame-n_fore+1:idx_frame
        % Ensure valid index
        if it_frame <= 0
            it_frame = 1;
        end
        for it_stream=1:2
            [header temp]=bkmload([folderName filesep streamName{select_img(it_stream)}],it_frame);
            temp_data{it_stream} = hilbert(double(temp(area_of_interest.axial_init:area_of_interest.axial_end,area_of_interest.lateral_init:area_of_interest.lateral_end)));
            temp_data{it_stream} = filtfilt(bpf.filter,temp_data{it_stream});
            new_line_count(1,it_stream)=header.lineCount;
        end
        PI_img_fore = abs(temp_data{select_img(1)}+temp_data{select_img(2)});
        PI_img_fore = interp2(X,Y,PI_img_fore,Xq,Yq,interpolation_type);
        new_time_stamp=header.timeStamp;
        
        fore_grnd_img = fore_grnd_img + PI_img_fore;
    end
    % Average by n_fore
    fore_grnd_img = fore_grnd_img/n_fore;
    
    % Remove background
    img = fore_grnd_img-bck_grnd_img;
    % Output img
    img = img.*(img >= 0);

elseif strcmp('B-mode',img_type) % For B-mode
    % Choose which images to use for B-mode
    select_img = [1];
    
    % Sum of n_fore images
    for it_frame = idx_frame-n_fore+1:idx_frame
        % Ensure valid index
        if it_frame <= 0
            it_frame = 1;
        end
        [header temp]=bkmload([folderName filesep streamName{select_img}],it_frame);
        temp_data = hilbert(double(temp(area_of_interest.axial_init:area_of_interest.axial_end,area_of_interest.lateral_init:area_of_interest.lateral_end)));
        temp_data = filtfilt(bpf.filter,temp_data);
        new_line_count = header.lineCount;

        B_img_fore = abs(temp_data{select_img(1)});
        B_img_fore = interp2(X,Y,B_img_fore,Xq,Yq,interpolation_type);
        new_time_stamp=header.timeStamp;
        
        fore_grnd_img = fore_grnd_img + B_img_fore;
    end
    % Average by n_fore
    fore_grnd_img = fore_grnd_img/n_fore;
    
    % Output img
    img = fore_grnd_img;

else
    erf('Not valid image type');
end


