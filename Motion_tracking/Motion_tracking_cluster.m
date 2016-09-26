%% Movement tracking B-mode
% B_mode
img_wind_cord = zeros(1,4);
idx_wind = 1;

img_wind_cord(idx_wind,:) = [350   550   114   116]


% Max movement y,x
max_mov_y = 10;
max_mov_x = 3;

% Setup of block_matching
fps = 50;
frames = 250;
start_frame = 1;
H = block_matching;
H.max_mov_y = max_mov_y;
H.max_mov_x = max_mov_x;
H.cost_function = 'xCorr';
H.img_wind_cord = img_wind_cord;
% mode = 'b_mode';

% Display Windows
% Dim and corr for boxes
for idx_wind = 1:size(img_wind_cord,1)               
    rectangle_corr(idx_wind,:) = [img_wind_cord(idx_wind,3), img_wind_cord(idx_wind,1), img_wind_cord(idx_wind,4)-img_wind_cord(idx_wind,3), img_wind_cord(idx_wind,2)-img_wind_cord(idx_wind,1)];
end

% Full image
img_disp = load_img_B_mode(1000);
figure(5); clf;
norm = max(abs(img_disp(:)));
limg=20*log10(abs(img_disp)/norm);
imagesc(limg,[-40 0]);% xlim([1 size(img,2)]); ylim([1 size(img,1)]);
colormap('gray'); xlabel('Lateral (mm)'); ylabel('Axial (mm)'); %title('B-mode image');
set(gca, 'DataAspectRatio',[1 1.67 1]) % set data aspect ratio in zoom box
set(gca, 'PlotBoxAspectRatio',[1 1 1])
%set(gca,'Xtick',linspace(0,280,5)); set(gca, 'XTickLabel',linspace(0,12,5));
%set(gca,'Ytick',linspace(0,1960,6)); set(gca, 'YTickLabel',linspace(0,25,6));

hold on;
for idx_wind = 1:size(img_wind_cord,1)
    rectangle('position',rectangle_corr(idx_wind,:),'EdgeColor','r','LineWidth', 2);
%     text(rectangle_corr(idx_wind,1),rectangle_corr(idx_wind,2)-50, int2str(idx_wind),'Color','r','FontSize',15,'FontWeight','bold');
end
set(gcf,'position',[-1850 570 560 420]);


%% Displacement             
vel_y = [];
vel_x = [];

% Run through all templates
for idx_wind = 1:size(img_wind_cord,1)
        % Load ref img
        idx_frame = start_frame;
        img_ref_wind = H.img_reference_window(idx_frame, idx_wind);
        % Repeat for all frames
        for idx_frame = start_frame:start_frame+frames-1;
            % Load ref img
            img_ref_wind = H.img_reference_window(idx_frame, idx_wind);
            
            % Load new img
            img_new_temp = H.img_new_template(idx_frame+1, idx_wind);%--
            
            [motion_y motion_x] = H.motion_displacement(img_ref_wind,img_new_temp);
            vel_y(idx_wind,idx_frame-start_frame+1) = motion_y;
            vel_x(idx_wind,idx_frame-start_frame+1) = motion_x;
        end        
end

%% Cluster part
% $Revision: 1.1 $  $Date: 02/13/14 13:10:27 $

% This is an example usage of running the tasks over the cluster.

% Define the directory where all the intermediate job handling files will
% be stored there
DataLocation = '/data/cfudata6/s134082/cluster_temp';
% cfu_mkdir(DataLocation);

% 'fcfu1' to 'fcfu12' or 'local'
sched = cluster_init('fcfu11',DataLocation);
% sched.ClusterMatlabRoot = '/usr/local/matlabR2010b_64bit';
global additionalSubmitArgs;
additionalSubmitArgs = '--cpus-per-task=2 --mem-per-cpu=400 --exclude=fcfu1,fcfu9 --nice=-10000';

% There is no limit on only assigning one job!
% you can assign as many jobs as you want and therefore assign many tasks
% to each of those jobs
job = createJob(sched);

% The workers on the cluster should know where to find the specified
% function, therefore one has to put the address of the directory in
% his/her startup.m file, which is normally at:
% /home/username/Documents/MATLAB/startup.m



%%
[X_comp,Y_comp] = meshgrid(15:20:175,150:100:850);


wind_size_y = 200;
wind_size_x = 3;
img_wind_cord = zeros(length(X_comp(:)),4);
for i = 1:length(X_comp(:))
    img_wind_cord(i,:) = [Y_comp(i)-wind_size_y/2,Y_comp(i)+wind_size_y/2,X_comp(i)-floor(wind_size_x/2),X_comp(i)+floor(wind_size_x/2)];
end
tasks = size(img_wind_cord,1);
for i=1:tasks
    disp(sprintf('Task:%d',i))
    createTask(job, @movement_tracking_b_mode_func, 3, {img_wind_cord(i,:), i});
end

%%
submit(job);
% Wait until all the tasks are finished then proceed to the next line

%%
wait(job);

%%
% Now that the job has finished we can harvest the results from all the
% workers, note that the results of the tasks are stored as a class object!
% data = fetchOutputs(job);
data = [];
for i=1:size(job.Tasks.OutputArguments,2)
    data{i}=job.Tasks(i,1).OutputArguments;
end
%%
% Clean the directory, this can be done at the end of your code whenever
% you don't need to assign anymore tasks to that specific job. Otherwise
% you still can assign new tasks to this job and run them on the cluster
destroy(job);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Job has finished and all the results have been read, the rest will be %
%  calculated on your local machine                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Check correct index
for i = 1:size(data,2)
    if isempty(data{1,i})
        disp(['Missing index ' num2str(i,'%d')])
%         data{1,i} = data{1,i-1};
%         data{1,i}{1,1}{3} = data{1,i}{1,1}{3}+1;
    end
%     if data{1,i}{1,1}{3} ~= i
%         disp(['Index ' num2str(i,'%d') ' not correct'])
%     end
end

vel_y_raw = [];
vel_x_raw = [];
for i = 1:size(data,2)
    vel_y_raw(i,:) = data{1,i}{1,1}{1};
    vel_x_raw(i,:) = data{1,i}{1,1}{2};
end

start_frame = 1;


