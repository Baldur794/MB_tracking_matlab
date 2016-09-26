function [ fore_grnd_img ] = load_img_contrast( idx_frame, idx_sync, n_bck, v_MB, mov_x_comp_contrast, mov_y_comp_contrast)
% Outputs foreground img.
filename = '/data/cfudata6/s134082/micro_bubble_data/mat_files/2016_04_21_15_29_41/frame_';

% img size from beamforming
% img_size = [684,123];% Contrast
img_size = [1, 1; 489, 61];% Contrast half removed
%img_size = [170, 5; 210, 28];% Small area with horizontal MB between 3970-4100
% img_size = [275, 15; 305, 23];% Small area with vertical MB between 4045-4100


interpolation_type = 'spline';
interpolation_factor_y = 5.1;
interpolation_factor_x = 19.7;
n_bck_grnd_skip = 50;
n_bck_jump = 10;
% contains foreground img and background img
bck_grnd_img = zeros(img_size(2,1)-img_size(1,1)+1, img_size(2,2)-img_size(1,2)+1);

if n_bck > 0
    % sum of n_bck_grnd images
    for i = idx_frame-n_bck_grnd_skip-(n_bck-1)*n_bck_jump:n_bck_jump:idx_frame-n_bck_grnd_skip
        load([filename num2str(i,'%d') '.mat'],'img');
        bck_grnd_img = bck_grnd_img + img(img_size(1,1):img_size(2,1),img_size(1,2):img_size(2,2));
    end
    % background average
    bck_grnd_img = bck_grnd_img/n_bck;
end

% Find normalized cut frequency
fps = 50;
MB_sigma = 0.2;
MB_size = 200*10^(-6);
x = -0.5:1/(MB_size/(v_MB*1/fps)):0.5;
gauss=1/sqrt(2*pi)/MB_sigma*exp(-(x).^2/2/MB_sigma/MB_sigma);
gauss_zero = [gauss zeros(1,1000-size(x,2))];
gauss_fft = fft(gauss_zero);
gauss_fft_norm = gauss_fft/max(abs(gauss_fft));
mag_cut = 0.1;
Wn_cut = find(abs(gauss_fft_norm) < mag_cut,1);
Wn_cut = Wn_cut/(size(gauss_fft_norm,2)/2);

% figure(30);plot(gauss);
% figure(31);plot(abs(gauss_fft_norm));

% Create filter
n_filter_1 = 100;
Wn_design = [0 Wn_cut-0.1*Wn_cut Wn_cut+0.1*Wn_cut 1];
mag = [1 1 0 0];
h_bandpass = fir2(n_filter_1,Wn_design,mag);
h_bandpass_flip = fliplr(h_bandpass);

% %----
% n_filter_2 = 20;
% % b = fir1(n,[0.05 [0.2 0.9]],'DC-1');
% b = fir1(n_filter_2,[0.05]);
% %----


[X,Y] = meshgrid(1:img_size(2,2)-img_size(1,2)+1,1:img_size(2,1)-img_size(1,1)+1);
[Xq,Yq] = meshgrid(1:1/interpolation_factor_x:img_size(2,2)-img_size(1,2)+1,1:1/interpolation_factor_y:img_size(2,1)-img_size(1,1)+1);

filt_img = zeros(size(Xq));
idx_filt = 1;
for idx_load = idx_frame-n_filter_1/2:idx_frame+n_filter_1/2
    idx_comp = idx_load-idx_sync+1;
    load([filename num2str(idx_load,'%d') '.mat'],'img');
    
    img = img(img_size(1,1):img_size(2,1),img_size(1,2):img_size(2,2))-bck_grnd_img;
    img = interp2(X,Y,img,Xq,Yq,interpolation_type);
%     figure(10123); imagesc(bck_grnd_img); colormap('gray');

    img_comp = img;
    
%     % compensated in axial direction
%     img_comp = zeros(size(img));
%     if mov_y_comp_contrast(idx_comp) < 0
%         img_comp(1:end+mov_y_comp_contrast(idx_comp),:)=img(-mov_y_comp_contrast(idx_comp)+1:end,:);
%     else
%         img_comp(mov_y_comp_contrast(idx_comp)+1:end,:)=img(1:end-mov_y_comp_contrast(idx_comp),:);
%     end
%     % compensated in lateral direction
%     if mov_x_comp_contrast(idx_comp) < 0
%         img_comp(:,1:end+mov_x_comp_contrast(idx_comp))=img_comp(:,-mov_x_comp_contrast(idx_comp)+1:end);
%     else
%         img_comp(:,mov_x_comp_contrast(idx_comp)+1:end)=img_comp(:,1:end-mov_x_comp_contrast(idx_comp));
%     end
%     
%     %----
%     img_comp = filter(b,1,img_comp,[],1);
%     img_comp = img_comp(n_filter_2/2+1:end,:);
%     %---

    filt_img = filt_img + img_comp*h_bandpass_flip(idx_filt);
    idx_filt = idx_filt + 1;
end
% fore_grnd_img = abs(filt_img-bck_grnd_img);

fore_grnd_img = abs(filt_img);
end

