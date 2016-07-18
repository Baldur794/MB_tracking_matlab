%% Load data
img = [];

n_start = 1;
n_end = 2;

for i = n_start:n_end
    img(:,:,i-n_start+1) = (load_img_B_mode(i));
end

%% Init block_matching
H = block_matching;


%%
[X Y] = H.center_position(img(:,:,1))

H.motion(real(img(:,:,1)),real(img(:,:,1)))
