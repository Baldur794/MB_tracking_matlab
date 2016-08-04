% folder='/data/cfudata7/cavh/phantoms/microbubble_BK/first_experiment/data/2016_04_21_14_32_21'; % Contrast

folder='/data/cfudata7/cavh/phantoms/microbubble_BK/first_experiment/data/2016_04_21_14_32_21'; % B-mode


% Contrast mode
% no_lines=123;
% file_size=168264;

% B-mode
 no_lines=560;
 file_size=2195200;

frame=0;
frame_skip=0;
total_frame_skip=0;
frame_skip_limit=10;
first=1;
run=true;

h_f=figure(1); clf

while run
    fname=sprintf('%s/frame_%d.bin', folder, frame);
    if ~exist(fname, 'file')
        frame_skip=frame_skip+1;
        total_frame_skip=total_frame_skip+1;
        if frame_skip==frame_skip_limit
            run=false;
        end
        continue;
    else
        frame_skip=0;
    end
    
    fprintf('Opening %s\n', fname);
    
    img=load_bk2300_zmq_data(fname, no_lines,(file_size/2)/no_lines, 'int16', 0);
    img = fliplr(img);
    env=abs(hilbert(img));
    if first
        norm=max(env(:));
        limg=20*log10(env/norm);
        drange=max(limg(:))+[-60 0];
        h_i=imagesc(limg, drange);
        h_t=title('Image');
        colormap gray
        first=0;
    else
        limg=20*log10(env/norm);
    end
    
    set(h_i, 'CData', limg);
    set(h_t, 'String', sprintf('Frame %d', frame));
    drawnow;
    
    frame=frame+1;
    pause(1/30)
    save(['/home/s134082/Desktop/Test/frame_' sprintf('%d',frame) '.mat'], 'img')
end

fprintf('Processed %d frames\n', frame-frame_skip_limit);
