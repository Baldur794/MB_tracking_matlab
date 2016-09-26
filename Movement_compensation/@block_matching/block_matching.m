classdef block_matching
    %BLOCK_MATCHING Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        max_mov_x = 25;
        max_mov_y = 25;
        cost_function = 'Norm Corr'
        img_wind_cord =[250,350,125,155];
        n_fore_grnd = 1;
        n_bck_grnd = 0;
        n_bck_grnd_skip = 1;
    end 
    
    methods(Access = public)
        % class constructor
        function [motion_y motion_x] = motion_displacement(obj, img_ref,img_new)
            if strcmp(obj.cost_function, 'Corr')
                % Matlab 2d NCC
                c = abs(xcorr2(img_ref*10^24,img_new*10^24));
                
                % Finds max correlation
                [max_y max_x] = peak_estimation(c, 'max');
                
                %                 % Finds max correlation
                %                 [max_y max_x] = find(c == max(c(:)));
                motion_y = max_y(1)-size(img_new,1)-obj.max_mov_y;
                motion_x = max_x(1)-size(img_new,2)-obj.max_mov_x;
                
            elseif strcmp(obj.cost_function, 'Norm Corr')
                % Matlab 2d NCC
                c = (normxcorr2(img_new*10^24,img_ref*10^24));
                c = c(size(img_new,1):size(c,1)-size(img_new,1)+1,size(img_new,2):size(c,2)-size(img_new,2)+1);
                
                % Finds peak position
                [max_y max_x] = peak_estimation(c, 'max');
                % Finds max correlation
                motion_y = max_y(1)-obj.max_mov_y-1;
                motion_x = max_x(1)-obj.max_mov_x-1;
            elseif strcmp(obj.cost_function, 'xCorr')
                for i = 1:obj.max_mov_x*2+1
                    c = [];
                    lag = [];
                    for j = 1:size(img_new,2)
                        [c(:,j) lag(:,j)] = crosscorr(img_new(:,j),img_ref(:,j+i-1),obj.max_mov_y*2);
                    end
                    c = c((size(c,1)-1)/2+1:(size(c,1)-1)/2+1+2*obj.max_mov_y,:);
                    c = mean(c',1)';
                    c_mean(:,i) = c;
                end
                [max_y max_x] = peak_estimation(c_mean, 'max');
                motion_x = max_x(1)-obj.max_mov_x-1;
                motion_y = max_y(1)-obj.max_mov_y-1;
%                 for i = 1:size(img_new,2)
%                     c = [];
%                     lag = [];
%                     for j = 1:obj.max_mov_x*2+1
%                         [c(:,j) lag(:,j)] = crosscorr(img_new(:,i),img_ref(:,j+i-1),obj.max_mov_y*2);
%                         %                      figure();crosscorr(img_new,img_ref(:,i),obj.max_mov_y*2);
%                     end
%                     c = c((size(c,1)-1)/2+1:(size(c,1)-1)/2+1+2*obj.max_mov_y,:);
%                     [max_y max_x] = peak_estimation(c, 'max');
%                     motion_x_temp(i) = max_x(1)-obj.max_mov_x-1;
%                 end
%                 motion_x = mean(motion_x_temp);
%                 
%                 c = [];
%                 lag = [];
%                 for i = 1:size(img_new,2)                 
%                     [c(:,i) lag(:,i)] = crosscorr(img_new(:,i),img_ref(:,obj.max_mov_x+i),obj.max_mov_y*2);
% %                     figure();crosscorr(img_new,img_ref(:,i),obj.max_mov_y*2);
%                 end
%                 c = c((size(c,1)-1)/2+1:(size(c,1)-1)/2+1+2*obj.max_mov_y,:);
%                 c = mean(c',1)';    
%                 [max_y max_x] = peak_estimation(c, 'max');
%                 motion_y = max_y(1)-obj.max_mov_y-1;
%                
%                 %----
% %                 c_axial_mean = mean(c',1);
% %                 [max_y max_x] = peak_estimation(c_axial_mean', 'max');
% %                 motion_y = max_y(1)-obj.max_mov_y-1;
%                 %----
                
            else
                % Cost matrix
                cost_matrix = zeros(size(img_ref,1)-size(img_new,1)+1,size(img_ref,2)-size(img_new,2)+1);
                
                for i = 1:size(cost_matrix,1)
                    for j = 1:size(cost_matrix,2)
                        y_start = i;
                        y_end = i+size(img_new,1)-1;
                        x_start = j;
                        x_end = j+size(img_new,2)-1;
                        img_ref_temp = img_ref(y_start:y_end,x_start:x_end);
                        if strcmp(obj.cost_function, 'MAD')
                            cost_matrix(i,j) = mad(img_ref_temp, img_new);
                        elseif strcmp(obj.cost_function, 'MSD')
                            cost_matrix(i,j) = msd(img_ref_temp, img_new);
                        else
                            error('Fail');
                        end
                    end
                end
                % Finds peak position
                [min_y min_x] = peak_estimation(cost_matrix,'min');
                
%                 [min_y min_x] = find(cost_matrix == min(cost_matrix(:)));
                motion_y = min_y-obj.max_mov_y-1;
                motion_x = min_x-obj.max_mov_x-1;
            end
        end
        
        function [img_ref_temp] = img_reference_window(varargin)   
            
             obj = varargin{1};
             idx_frame = varargin{2};
             idx_wind = varargin{3};
             if size(varargin,2) == 4
             mode = varargin{4};
             else
                 img_avg = varargin{4};
                 mode = varargin{5};
             end
             
            % Load img
            if strcmp(mode,'b_mode')
                img_ref = abs(load_img_B_mode(idx_frame));%---------------------------
            elseif strcmp(mode,'contrast')
                img_ref = abs(load_img_contrast(idx_frame,obj.n_fore_grnd,obj.n_fore_grnd,obj.n_bck_grnd_skip));
                img_ref(img_avg == 0) = 0;
            end
            % Find coordinates
            y_start = obj.img_wind_cord(idx_wind,1)-obj.max_mov_y;
            y_end = obj.img_wind_cord(idx_wind,2)+obj.max_mov_y;
            x_start = obj.img_wind_cord(idx_wind,3)-obj.max_mov_x;
            x_end = obj.img_wind_cord(idx_wind,4)+obj.max_mov_x;
            
            % Make template window
            img_ref_temp = zeros(y_end-y_start+1,x_end-x_start+1);
            img_ref_temp = img_ref(y_start:y_end,x_start:x_end);
        end
        
        function [img_new_temp] = img_new_template(varargin)
            obj = varargin{1};
            idx_frame = varargin{2};
            idx_wind = varargin{3};
            if size(varargin,2) == 4
                mode = varargin{4};
            else
                img_avg = varargin{4};
                mode = varargin{5};
            end
            
            % Load new img
            if strcmp(mode,'b_mode')
                img_new = abs(load_img_B_mode(idx_frame));%-------------------
            elseif strcmp(mode,'contrast')
                img_new = abs(load_img_contrast(idx_frame,obj.n_fore_grnd,obj.n_fore_grnd,obj.n_bck_grnd_skip));
                img_new(img_avg == 0) = 0;
                %abs(load_img_contrast(idx_frame,obj.n_fore_grnd,obj.n_fore_grnd,obj.n_bck_grnd_skip));
            end
            
            % Find coordinates
            y_start = obj.img_wind_cord(idx_wind,1);
            y_end = obj.img_wind_cord(idx_wind,2);
            x_start = obj.img_wind_cord(idx_wind,3);
            x_end = obj.img_wind_cord(idx_wind,4);
            
            img_new_temp = zeros(y_end-y_start+1,x_end-x_start+1);          
            img_new_temp = img_new(y_start:y_end,x_start:x_end);
        end
    end
end

%% Subfunctions

function mad_cost = mad(img_ref_temp, img_new)
    diff = abs(img_ref_temp-img_new);
    mad_cost = sum(diff(:))/length(img_new(:));
end

function msd_cost = msd(img_ref_temp, img_new)
    diff = (img_ref_temp-img_new).^2;
    msd_cost = sum(diff(:))/length(img_new(:));
end

function [peak_y peak_x] = peak_estimation(match_matrix, min_max)
    if strcmp(min_max, 'max')
        [y_m x_m] = find(match_matrix == max(match_matrix(:)));
    elseif strcmp(min_max, 'min')
        [y_m x_m] = find(match_matrix == min(match_matrix(:)));
    else
        error('Error')
    end
    x_m = x_m(1);
    y_m = y_m(1);
    if (x_m == size(match_matrix,2)) || (x_m == 1)
        peak_x = x_m;
    else
        peak_x = x_m-(match_matrix(y_m,x_m+1)-match_matrix(y_m,x_m-1))/(2*(match_matrix(y_m,x_m+1)-2*match_matrix(y_m,x_m)+match_matrix(y_m,x_m-1)));
    end
    if (y_m == size(match_matrix,1)) || (y_m == 1)
        peak_y = y_m;
    else
        peak_y = y_m-(match_matrix(y_m+1,x_m)-match_matrix(y_m-1,x_m))/(2*(match_matrix(y_m+1,x_m)-2*match_matrix(y_m,x_m)+match_matrix(y_m-1,x_m)));
    end
end





