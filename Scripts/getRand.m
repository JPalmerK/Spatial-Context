function AdjRand= getRand(obj)

% Check if cluster id's are present if not update it
if isempty(obj.Cluster_id)
    updateClusterID(obj)
end


% Pull out the arrival array for easy reading
arrival_array = gather(obj.arrivalArray);
trueClusters =unique(arrival_array(:,end));


% Create acoustic encounters
AcousticEncounters = acEnc(obj);
AcousticEncIds = unique(AcousticEncounters);

% True agents
Truth = (arrival_array(:,end));

NewTruth =(ones(size(Truth)));

TruthClusterId =1;

nextIdx = max(AcousticEncIds)+1;
for ii=1:length(AcousticEncIds)

    % Index of the acoustic encounter
    encounter_idx = find(AcousticEncounters==AcousticEncIds(ii));
    
    % Agent ID's within that acoustic encounter
    agentIDs = Truth(encounter_idx);
    
    uniqueAgents = unique(agentIDs);
    newClusterIDs = ones(size(encounter_idx));
    for jj=1:length(uniqueAgents)
        
        newClusterIDs(agentIDs == uniqueAgents(jj)) =  TruthClusterId;
        TruthClusterId=TruthClusterId+1;
        
    end
        % Update the truth clusters
        NewTruth(encounter_idx) = newClusterIDs;
        
    
end


    % Get the adjusted rand (third party)
    [f,~,~,~] = RandIndex(obj.Cluster_id, NewTruth);
    AdjRand = f;


end