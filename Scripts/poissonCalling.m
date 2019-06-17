%% Bout information
% Details from Mathews et al 2001
% intervals between call clusters 0.57-1.63 min
ave_z = [1.51, 0.57,0.21, 0.6, 0.63];

% Spacing between each call in each cluster
ave_call_spacing = [0.33,0.33,0.17,0.25,0.2];

% Number of calls in each cluster
ave_size = median([1.03,1.43,1.54,1.96,4.13]);

% Calls produced in bouts, nested poisson distribution not
% implemented

sim_duration = 60*60


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1) Pick the number of clusters

% pick a random index from the reported values
id = randi(length(ave_size));

% Rate of clusters per minute
cluster_rate = 1/ave_z(id);
n_clusters = poissrnd(sim_duration/60/cluster_rate);


% 2) Pick the number of calls in each cluster
cluster_size = poissrnd(max(ave_size), [1 n_clusters]);

% Need at least one call in the cluster
cluster_size = max([ones(size(cluster_size)); cluster_size]);


% 3) Determine where in the record the calls are

% create start times for all calls
call_produced = zeros(1, sum(cluster_size));
call_id =1;


% Iterate through the clusters and get the start time of each cluster and
% subsiquently each call.

for ii =1:n_clusters
    
    % Create an average spacing between the clusters
    cluster_space = gamrnd(1, 181, 1);
    
    % space between the calls in seconds
    call_space = [0 gamrnd(1, 5, ...
        [1, cluster_size(ii)-1])];
    
    % Start of each call sum of the call spacings plus the start of the
    % cluster plus the end of the last call
    call_produced(call_id:(call_id+ cluster_size(ii)-1)) = ...
        cumsum(call_space)+cluster_space + max(call_produced);
    
    % iterator
    call_id= call_id+length(call_space);
    
    
end

% Remove any calls with start times after the end of the simulation
% duration

call_produced = sort(call_produced(call_produced<sim_duration));


% Plotting
figure(1)
subplot(3,1,1)
scatter(1:length(call_produced), call_produced, '.')
subplot(3,1,2)
scatter(call_produced, ones(size(call_produced)), '.')

subplot(3,1,3)
hist(diff(call_produced),50)



