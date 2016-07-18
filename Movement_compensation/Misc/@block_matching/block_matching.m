classdef block_matching
    %BLOCK_MATCHING Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        block_size_x = 11;
        block_size_y = 11;
        max_disp_x = 10;
        max_disp_y = 10;
        overlap_x = 0;
        overlap_y = 0;
    end 
    
    methods(Access = public)
    % class constructor
        function [center_pos_x center_pos_y] = center_position(obj,img)
           [center_pos_x center_pos_y] = center_position(obj,img)
        end
        
        function  motion(obj,img_ref,img)
           motion_displacement(obj,img_ref,img)
        end
    end
end

%% Subfunctions
function [center_pos_x center_pos_y] = center_position(obj,img)
    center_pos_x = zeros(1+round((size(img,1)-obj.block_size_y)/(obj.block_size_y-obj.overlap_y)),...
        1+round((size(img,2)-obj.block_size_x)/(obj.block_size_x-obj.overlap_x)));
    center_pos_y = zeros(size(center_pos_x));
    
    for i = 1:size(center_pos_x,2)
        center_pos_x(:,i) = ceil(obj.block_size_x/2) + (obj.block_size_x - obj.overlap_x) * (i-1);
    end
    
    for i = 1:size(center_pos_y,1)
        center_pos_y(i,:) = ceil(obj.block_size_y/2) + (obj.block_size_y - obj.overlap_y) * (i-1);
    end
    
end


function motion = motion_displacement(obj,img_ref,img)
    % Init windows
    img_ref_temp = zeros(obj.block_size_y,obj.block_size_x);
    img_temp = zeros(obj.block_size_y,obj.block_size_x);
    
    % Matrix with center pos
    [center_pos_x center_pos_y] = center_position(obj,img);
    
    % Vectors with center pos
    center_x = center_pos_x(1,:);
    center_y = center_pos_y(:,1)';
    
    for i = 1:length(center_y)
        for j = 1:length(center_x)
            
            y_start = center_y(i) - floor(obj.block_size_y/2);
            y_end = center_y(i) + floor(obj.block_size_y/2);
            x_start = center_x(j) - floor(obj.block_size_x/2);
            x_end = center_x(j) + floor(obj.block_size_x/2);
            
            img_ref_temp = img(y_start:y_end,x_start:x_end);         
            
            corr = zeros(obj.max_disp_y*2+1, obj.max_disp_x*2+1);
            for k = -obj.max_disp_y:obj.max_disp_y
                for h = -obj.max_disp_x:obj.max_disp_x
                    
                    y_start_temp = y_start + k;
                    y_end_temp = y_end + k;
                    x_start_temp = x_start + h;
                    x_end_temp = x_end + h;
                    
                    if ( y_start_temp < 1 || y_end_temp > size(img,1) ...
                        || x_start_temp < 1 || x_end_temp > size(img,2))
                        continue;
                    end
                    
                    img_temp = img(y_start_temp:y_end_temp,x_start_temp:x_end_temp);
                    
                    c = normxcorr2(img_ref_temp,img_temp);
                    corr(k+obj.max_disp_y+1,h+obj.max_disp_x+1) = c(obj.max_disp_y+1,obj.max_disp_x+1);
                               
                end
            end
            
        end
    end
end


% c = normxcorr2(img_template,img_new);
% corr_fix(idx_temp,idx_seq) = c(img_template_corr(idx_temp,2),img_template_corr(idx_temp,4));
% [max_y max_x] = find(c == max(c(:)));
% corr_max_y(idx_temp,idx_seq) = max_y;
% corr_max_x(idx_temp,idx_seq) = max_x;