%% Acoustic Encounter -
%Cluster based only on time of arrivals Baseline (step 4) equivallent to acoustic encounters
function cluster_vals = acEnc(obj)


% Look for gaps bigger than the maximum overlap time
diff_vals = diff(obj.arrivalArray(:,1));

%             % Find gaps bigger than the allowable time
%             time_idx = quantile(diff_vals, .95);
%

%%%%%%%%%%%%%%%%%%%%%%%%
% For sensitivity analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_idx = obj.maxEltTime;
% Find all break points
idx = find(diff_vals>time_idx);

% Double check!! 
cluster_vals= ones(size(obj.arrivalArray(:,1)));
cluster_id =2;


for ii = 1:length(idx)
    
    cluster_vals(idx(ii)+1:end) = cluster_id;
    cluster_id = cluster_id +1;
    
end



end