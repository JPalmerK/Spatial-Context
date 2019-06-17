function newClusterId = allignclusters(agent_ids, Cluster_id)
% Function to allign clusters
% agent_ids - true agent clusters from simulation
% Cluster_id - cluster ids created by the clustering algorithim



refed_cluster_id = zeros(size(Cluster_id))/0;

group_start = 1;




% Step through the known clusters and add the reference
for ii = 1:length(refed_cluster_id)
    
    
    if isnan(refed_cluster_id(ii)) && ...
            group_start<length(unique(agent_ids));
        
        linked_cluster_id = Cluster_id(ii);
        refed_cluster_id(find(Cluster_id == linked_cluster_id)) = group_start;
        group_start =group_start+1;
        
        
    end
    
    
    
end


% Get the remaining Nan Values from the refed_cluster_id

% Step through the known clusters and add the reference
for ii = 1:length(Cluster_id)
    
    if isnan(refed_cluster_id(ii))
        
        linked_cluster_id = Cluster_id(ii);
        refed_cluster_id(find(Cluster_id == linked_cluster_id)) = group_start;
        group_start =group_start+1;
        
        
    end
    
    
    
end


newClusterId = refed_cluster_id;


end