function [ fore_grnd_img ] = load_img_contrast( idx_seq, n_bck_grnd, n_bck_grnd_skip, v_MB)
% Outputs foreground img.
filename = '/data/cfudata6/s134082/Bachelorprojekt/micro_bubble_data/mat_files/2016_04_21_15_29_41/frame_';

% img size from beamforming
%img_size = [684,123];% Contrast
img_size = [489,61];% Contrast half removed

% contains foreground img and background img
bck_grnd_img = zeros(img_size);

if n_bck_grnd > 0
    % sum of n_bck_grnd images
    for i = idx_seq-n_bck_grnd_skip*n_bck_grnd-n_bck_grnd_skip:n_bck_grnd_skip:idx_seq-2*n_bck_grnd_skip
        load([filename num2str(i,'%d') '.mat'],'img');
        bck_grnd_img = bck_grnd_img + img(1:img_size(1),1:img_size(2));
    end
    % background average
    bck_grnd_img = bck_grnd_img/n_bck_grnd;
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
Wn_cut = find(abs(gauss_fft_norm)<mag_cut,1);
Wn_cut = Wn_cut/(size(gauss_fft_norm,2)/2);

% figure(30);plot(gauss);
% figure(31);plot(abs(gauss_fft_norm));

% Create filter
n_filter = 100;
Wn_design = [0 .05 .055 Wn_cut-0.1*Wn_cut Wn_cut+0.1*Wn_cut 1];
mag = [1 1 1 1 0 0];
h_bandpass = fir2(n_filter,Wn_design,mag);
h_bandpass_flip = fliplr(h_bandpass);

filt_img = zeros(img_size);
idx_filt = 1;
for i = idx_seq-n_filter/2:idx_seq+n_filter/2
    load([filename num2str(i,'%d') '.mat'],'img');
    filt_img = filt_img + img(1:img_size(1),1:img_size(2))*h_bandpass_flip(idx_filt);
    idx_filt = idx_filt + 1;
end
fore_grnd_img = abs(filt_img-bck_grnd_img);
%fore_grnd_img = abs(filt_img);
end

