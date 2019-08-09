function truthAndpreds =createSpeciesPreds(obj)
% This function is the classifier function that could use a bit
% of love. It's a bit hacky. Looks at calls and assigns a
% random probability score indicating that each call came from
% species 1.

% Check the arrival array, if it's empty populate it
if isempty(obj.arrivalArray)
    UpdateArrArray(obj);
end
% Check to see if clusters id had been added
if isempty(obj.Cluster_id)
    updateClusterID(obj);
end



% Stick the real id's next to the predicted clusters
truthAndpreds =[obj.arrivalArray(:,end) obj.Cluster_id];
truthAndpreds(:,3)=0; % True class
truthAndpreds(:,4)=0; % Score
truthAndpreds(:,5)=0; % Prediction after applying likelihood ratio

truthAndpreds= array2table(truthAndpreds, 'VariableNames',...
    {'TrueClust','PredClust', 'TrueSpp', 'Score',...
    'ClusterScore'});


% Detector Parameters
rw_mean = obj.betaParm1;
rw_sd = obj.betaParm2;

% For each detected agent, assign a species and classification
% probability
agent_idxs =unique(truthAndpreds.TrueClust);

for ii =1:length(agent_idxs)
    
    % Agent id
    agent_id =agent_idxs(ii);
    
    % If it's odd it's definitely a humpback
    if mod(agent_id,2)
        
        % Get the indicies of the calls of that individual
        rw_idx = find(truthAndpreds.TrueClust == agent_id);
        
        % Pull from distribution to estimate classificaiton probability
        score = betarnd(obj.betaParm1, obj.betaParm2,[1 length(rw_idx)]);
        %score = .6;
        
        
        truthAndpreds.TrueSpp(rw_idx) =1; %yes target spp.!
        truthAndpreds.Score(rw_idx) = score';
        
        
    else
        
        
        % Same as above, only other species
        
        % Get the indicies of the calls of that individual
        mn_idx = find(truthAndpreds.TrueClust == agent_id);
        
        % Transform so that it's between 0 and 1
        score = 1- betarnd(obj.betaParm1, obj.betaParm2,[1 length(mn_idx)]); % logit (or inverse logit)
        %score = .4;
        
        % Update the species (binary) and the classifier
        % predication score
        truthAndpreds.TrueSpp(mn_idx) = 0; % No! Not target spp.
        truthAndpreds.Score(mn_idx) = score';
        
    end
    
    
    
end



end
