function cluster = updateChainsEncounterFirst(obj)
% Use the similarity matrix to update the cluster chains



% Ladder Linkages rebuild
simmat= obj.Sim_mat;
arraivalArray = gather(obj.arrivalArray);
max_gap = obj.maxEltTime;
simthresh = obj.cutoff;

% Prune the simulation matrix such that each row contains only members of
% the acoustic encounter

for ii= 1:size(simmat,1)

    % Get the diff times
    diff_times = diff(arraivalArray(ii:end,1));
    
    break_idx = find(diff_times>max_gap,1);
    
    simmat(ii, (ii+break_idx-1):end) = NaN;
     simmat((ii+break_idx-1):end, ii) = NaN;

end


clusterN=1;
cluster = struct();

% Index of where we are in the array
idx = 1:size(arraivalArray,1);

% For each row, find the first gap
while size(simmat,1)>1
    
    col_idx = idx(1);
    
    simscores = simmat(1,:);
    % If the the first gap larger than the maximum elapsed time it gets
    % it's own cluster
    
    % otherwise get the scores and remove the ones that are above the threshold
    temp=1:length(simscores);
    matching_idx = unique([1, temp(simscores>=simthresh)]);
    
    
    cluster(clusterN).index = col_idx+matching_idx-1;
    cluster(clusterN).n = length(matching_idx);
    
    
    
    % Remove the rows from thematrix
    simmat(matching_idx,:)=[];
    simmat(:,matching_idx)=[];
    idx(matching_idx)=[];
    
    
    clusterN=clusterN+1;
end



end





