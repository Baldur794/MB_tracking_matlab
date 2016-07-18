%% Load Img's
img = [];

n_start = 1;
n_end = 300;

for i = n_start:n_end
    img(:,:,i-n_start+1) = (load_img_B_mode(i));
end
%% Motion calculation
% norm = max(img_B_mode(:));
% limg=20*log10(img_B_mode/norm);
%figure(1);
%imagesc(limg,[-40 0]); colormap(gray);


block_size_x = 9;
block_size_y =19; 
max_disp_x = 10;
max_disp_y = 10;

H = vision.BlockMatcher('ReferenceFrameSource', 'Input port', 'BlockSize', [block_size_y block_size_x],...
    'MaximumDisplacement', [max_disp_y max_disp_x], 'SearchMethod', 'Exhaustive', 'Overlap', [0 0],'MatchCriteria','Mean square error (MSE)');
H.OutputValue = 'Horizontal and vertical components in complex form';


motion = [];
for i = 1:size(img,3)-1
    motion(:,:,i) = step(H, real(img(:,:,1)), real(img(:,:,i+1)));
end


 
 %% Motion video
 
[X, Y] = meshgrid(1:block_size_x:size(img(:,:,1), 2), 1:block_size_y:size(img(:,:,1), 1));
 
norm=max(img(:));

outputVideo=VideoWriter('Motion_movement.avi');
outputVideo.FrameRate=4;
open(outputVideo);

mov(1:size(motion,3))= struct('cdata',[],'colormap',[]);
figure(3)

norm = max(abs(img(:)));
limg=20*log10(abs(img(:,:,1))/norm);
imagesc(limg,[-40 0]);% xlim([1 size(img,2)]); ylim([1 size(img,1)]);
colormap('gray'); xlabel('Lateral (mm)'); ylabel('Axial (mm)'); %title('B-mode image');
set(gca, 'DataAspectRatio',[1 3.46 1]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])
set(gca,'Xtick',linspace(0,280,5)); set(gca, 'XTickLabel',linspace(0,12,5));
set(gca,'Ytick',linspace(0,1960,6)); set(gca, 'YTickLabel',linspace(0,25,6));
hold on;

pause(0.3);
for i = 1:size(motion,3)
    motion_temp = motion(:,:,i);
    
    limg=20*log10(abs(img(:,:,i))/norm);
    imagesc(limg,[-40 0]);% xlim([1 size(img,2)]); ylim([1 size(img,1)]);
    colormap('gray'); xlabel('Lateral (mm)'); ylabel('Axial (mm)'); %title('B-mode image');
    set(gca, 'DataAspectRatio',[1 3.46 1]) % set data aspect ratio in zoom box
    set(gca, 'PlotBoxAspectRatio',[1 1 1])
    set(gca,'Xtick',linspace(0,280,5)); set(gca, 'XTickLabel',linspace(0,12,5));
    set(gca,'Ytick',linspace(0,1960,6)); set(gca, 'YTickLabel',linspace(0,25,6));
    hold on;
    
    quiver(X(:), Y(:), real(motion_temp(:)), imag(motion_temp(:)), 0);
    mov=getframe(gcf);
    writeVideo(outputVideo,mov.cdata);
    pause(0.3);
end
close(gcf)
close(outputVideo);


%%
motion_point = squeeze((motion(41:43,17,:)));
figure();
motion_point_test = mean(motion_point);
plot(imag(motion_point)');
figure();
plot(imag(motion_point_test));
sum(motion_point_test)

%%

img_temp_comp=img(700:900,130:160,:);
img_temp=img(700:900,130:160,:);
for i = 1:size(img_temp_comp,3)-1
    img_temp_comp(:,:,i+1)= imtranslate(img_temp_comp(:,:,i+1), [0 imag(motion_point_test(i))]);
end
%%


outputVideo=VideoWriter('Motion_movement_temp.avi');
outputVideo.FrameRate=1;
open(outputVideo);

mov(1:size(motion,3))= struct('cdata',[],'colormap',[]);
figure(3)

norm = max(abs(img_temp_comp(:)));
pause(0.3)

for i = 1:50
    
    
    limg=20*log10(abs(img_temp(:,:,i))/norm);
    imagesc(limg,[-40 0]);% xlim([1 size(img,2)]); ylim([1 size(img,1)]);
    colormap('gray'); xlabel('Lateral (mm)'); ylabel('Axial (mm)'); %title('B-mode image');
    set(gca, 'DataAspectRatio',[1 3.46 1]) % set data aspect ratio in zoom box
    set(gca, 'PlotBoxAspectRatio',[1 1 1])
    %set(gca,'Xtick',linspace(0,280,5)); set(gca, 'XTickLabel',linspace(0,12,5));
    %set(gca,'Ytick',linspace(0,1960,6)); set(gca, 'YTickLabel',linspace(0,25,6));
    
    
    mov=getframe(gcf);
    writeVideo(outputVideo,mov.cdata);
    pause(0.3);
end
close(gcf)
close(outputVideo);
