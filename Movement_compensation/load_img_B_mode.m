function [ fore_grnd_img ] = load_img( idx_seq )
% Outputs foreground img.
filename = '/data/cfudata6/s134082/Bachelorprojekt/micro_bubble_data/mat_files/2016_04_21_14_32_21/frame_';
% filepath = '/data/cfudata6/s134082/Bachelorprojekt/simulation_data/spinning_disk/Spinningdisk_0.002m_s/FIELD/bf_em_data/flow_BF_res1010_full_new/seq_';

% chosen emission
em_idx = 1;

% img size from beamforming
%img_size = [1960,560];% B-mode 
img_size = [1960,280];% B-mode half removed 

% contains foreground img and background img
fore_grnd_img = zeros(img_size);
bck_grnd_img = zeros(img_size);

% % background options
% n_bck_grnd = 1000;% # of img's to make back ground
% n_bck_grnd_skip = 1;% # of img's skipped between each background img

% % sum of n_bck_grnd images
% for i = idx_seq-n_bck_grnd_skip*n_bck_grnd:n_bck_grnd_skip:idx_seq-n_bck_grnd_skip
%     %bck_grnd_img = bck_grnd_img + cell2mat(struct2cell(load([filename num2str(i,'%04d') '/em_' num2str(em_idx,'%04d') '.mat'],'bf_data')));
%     load([filename num2str(i,'%d') '.mat'],'img');
%     bck_grnd_img = bck_grnd_img + img;
% end
% % background average
% bck_grnd_img = bck_grnd_img/n_bck_grnd;

% foreground from present frame subtracted the background
%fore_grnd_img = abs(cell2mat(struct2cell(load([filename num2str(idx_seq,'%04d') '/em_' num2str(em_idx,'%04d') '.mat'],'bf_data'))) - bck_grnd_img);

load([filename num2str(idx_seq,'%d') '.mat']);
fore_grnd_img = (hilbert(img(1:2:end,1:img_size(2))));% - bck_grnd_img;

%  load([filepath num2str(idx_seq,'%04d') '/em_0003.mat']); 
%  fore_grnd_img = bf_data;
end

