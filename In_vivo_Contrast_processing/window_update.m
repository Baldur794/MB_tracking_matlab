function [ MB_window_coord MB_window_out_of_bounce] = window_update( MB, MB_index, MB_window_1, MB_window_2, img_size, age)
% Returns window coordinates and a flag indicating out of img.

% Update window coordinates
if age == 1
    MB_window_coord(1,1) = MB(MB_index).new_pos(1)-MB_window_1(1);
    MB_window_coord(2,1) = MB(MB_index).new_pos(1)+MB_window_1(1);
    MB_window_coord(1,2) = MB(MB_index).new_pos(2)-MB_window_1(2);
    MB_window_coord(2,2) = MB(MB_index).new_pos(2)+MB_window_1(2);
else
    MB_window_coord(1,1) = MB(MB_index).new_pos(1)-MB_window_2(1)+MB(MB_index).vel(1);
    MB_window_coord(2,1) = MB(MB_index).new_pos(1)+MB_window_2(1)+MB(MB_index).vel(1);
    MB_window_coord(1,2) = MB(MB_index).new_pos(2)-MB_window_2(2)+MB(MB_index).vel(2);
    MB_window_coord(2,2) = MB(MB_index).new_pos(2)+MB_window_2(2)+MB(MB_index).vel(2);
end

% Check for out of bounce
% y_start
MB_window_out_of_bounce = 0;
if MB_window_coord(1,1) <= 0
    MB_window_coord(1,1) = 1;
    MB_window_out_of_bounce = 1;
end
% y_end
if MB_window_coord(2,1) > img_size(1)
    MB_window_coord(2,1) = img_size(1);
    MB_window_out_of_bounce = 1;
end
% x_start
if MB_window_coord(1,2) <= 0
    MB_window_coord(1,2) = 1;
    MB_window_out_of_bounce = 1;
end
% x_end
if MB_window_coord(2,2) > img_size(2)
    MB_window_coord(2,2) = img_size(2);
    MB_window_out_of_bounce = 1;
end

end

