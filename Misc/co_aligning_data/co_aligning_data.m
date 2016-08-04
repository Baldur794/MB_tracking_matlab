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


n = (1:20);
test1 = zeros(1,20);%sin(2*pi/10*n)+sin(2*pi/5*n);
test1(10) = 1;

test2 = filter(h_bandpass,1,test1);

idx_seq = 4000;
test3 = [];
idx_filt = 1;
for i = 1:n_filter+1
%    load([filename num2str(i,'%d') '.mat'],'img');
%     filt_img = filt_img + img(1:img_size(1),1:img_size(2))*h_bandpass(idx_filt);
%     filt_img = filt_img + (img(1:img_size(1),1:img_size(2))-bck_grnd_img)*h_bandpass(idx_filt);
    test3 = test3 + test1(i)*h_bandpass(idx_filt); 
    idx_filt = idx_filt + 1;
end
%%
data = [4 1 2 2 1 3 5];
h11 = [0 2 0 1 0];
n_filter = size(h11,2)-1
test11111 = filter(h11,1,test11)
h111 = fliplr(h11);
idx = 4;
test1111 = 0;
idx_filter = 1;
for i = idx-n_filter/2:idx+n_filter/2
    test1111 = test1111 + h111(idx_filter)*test11(i);
    idx_filter = idx_filter +1;
end
test11111(n_filter/2+1:end)
test1111









