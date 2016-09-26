%% Cluster
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
additionalSubmitArgs = '--cpus-per-task=2 --mem-per-cpu=400 --exclude=fcfu1,fcfu5,fcfu6,fcfu9 --nice=-10000';

% There is no limit on only assigning one job!
% you can assign as many jobs as you want and therefore assign many tasks
% to each of those jobs
job = createJob(sched);

% The workers on the cluster should know where to find the specified
% function, therefore one has to put the address of the directory in
% his/her startup.m file, which is normally at:
% /home/username/Documents/MATLAB/startup.m

%%

idx_frame_start = 1000:1000:37000;

for i=1:size(idx_frame_start,2)
        disp(sprintf('Task:%d',i));
        createTask(job, @MB_tracking_func_4, 2, {idx_frame_start(i), i});
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
for i=1:size(job.Tasks,1)
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
%%
MB_log = [];
for i = 1:1
    load(['Task' num2str(i,'%d') '.out.mat'])
    MB_log = [MB_log argsout{1,1}{1,1}];
end


