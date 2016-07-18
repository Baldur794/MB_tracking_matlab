function [ mean_intensity ] = mean_intensity_3d( data_matrix_3d )
% Calculates mean intensity of 3D matrix
%   [x y z] = weighted_centroid_3d(data_maxtrix_3D)
mean_intensity = sum(data_matrix_3d(:))/length(data_matrix_3d(data_matrix_3d > 0));
end

