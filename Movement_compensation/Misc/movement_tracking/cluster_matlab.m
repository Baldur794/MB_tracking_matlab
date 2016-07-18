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
[X,Y] = meshgrid(25:40:185,150:700:1550);

max_mov_y = 16;
max_mov_x = 3;
wind_size_y = 200;
wind_size_x = 20;

img_wind_cord = zeros(length(X(:)),4);
for i = 1:length(X(:))
    img_wind_cord(i,:) =[Y(i)-wind_size_y/2,Y(i)+wind_size_y/2,X(i)-wind_size_x/2,X(i)+wind_size_x/2];
end
tasks = length(X(:));
for i=1:tasks
    disp(sprintf('Task:%d',i))
    createTask(job, @movement_tracking_b_mode_func, 3, {img_wind_cord(i,:), 1});
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
for i=1:tasks
    data{i}=job.Tasks(i,1).OutputArguments;
end
% Clean the directory, this can be done at the end of your code whenever
% you don't need to assign anymore tasks to that specific job. Otherwise
% you still can assign new tasks to this job and run them on the cluster
destroy(job);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Job has finished and all the results have been read, the rest will be %
%  calculated on your local machine                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
mov_y = [];
mov_x = [];
for i = 1:tasks
    mov_y(i,:) = data{1,i}{1,1}{1};
    mov_x(i,:) = data{1,i}{1,1}{2};
end