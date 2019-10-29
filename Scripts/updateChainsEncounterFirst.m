function chain = updateChainsEncounterFirst(obj)
% Use the similarity matrix to update the cluster chains


% Use the similarity matrix to update the cluster chains

% If similarity matrix was pre-allocated, recreate it

if ismatrix(obj.Sim_mat)
% Pull out a bunch of the variables to make life easier
simmat= (obj.Sim_mat); % The similarity matrix
else

    % Pull out the similarity values for each call
    
end




arrivalArray = (obj.arrivalArray); % Arrival times for the calls
max_gap = obj.maxEltTime; % user parameter
simthresh = obj.cutoff;% user parameter

% Initialzie some values
ids = 1:size(arrivalArray,1); % origional ID's fo the calls 
timeDiffs = [0; diff(arrivalArray(:,1))]; % Time diffs for acoustic encs.
% Indicies fo the start of each acoustic encounter
acousticEncounterBreaks = [1; (find(timeDiffs>max_gap))];

% Initialize a datatable for indexing around and preserving the origional
% information
dataTable = table();
dataTable.timeDiffs = timeDiffs;
dataTable.TrueIDs = ids';
dataTable.AcousticEncounters = ids'./ids';% this is just zero
dataTable.Clustered = ids'-ids'; % this is just zero
dataTable.ArrivalTimes = arrivalArray(:,1); % Arrival time on the parent
dataTable.TrueCluster = arrivalArray(:,end); % Dummy, unused for testing

% Demarkate the acoustic encounters in the datatable
for ii=1:length(acousticEncounterBreaks)
    
    dataTable.AcousticEncounters(acousticEncounterBreaks(ii):end) =ii;
    
end

% Initialzie the cluster id and chain structure 
clusterID =1;
chain =struct();


% current index into the encouter
idxInEnc =1;

% Step through acoustic encounters and define sub-clusters based on
% similarity scores and elapsed time
for ii=1:length(unique(dataTable.AcousticEncounters))
    
    % pull out the acoustic encounter
    encounterSub = dataTable(dataTable.AcousticEncounters==ii,:);
    encounterSub.SimScore=encounterSub.Clustered;
    encounterSub.ElapsedTime = encounterSub.ArrivalTimes- ...
        encounterSub.ArrivalTimes(idxInEnc);
    
    % Get the similartiy scores for first detection in the dataset
    encounterSub.SimScore(idxInEnc:end) = simmat(encounterSub.TrueIDs(idxInEnc),...
        encounterSub.TrueIDs(idxInEnc):...
        encounterSub.TrueIDs(end))';
    
    % Start of the acoustic encounter, the first call is always part of the
    % next cluster (initialize the 'Cluster' bool to 1 for the first
    % encounter
    encounterSub.Clustered(1)=1;
    
    % Index into the acoustic encounter
    idxInEnc =1;
    
    % Step throug the acoustic encounter making sub encounters where the
    % spatial contex matches up
    while height(encounterSub)>0
        
        % Identify all chains in the cluster- calls with elapsed times less
        % than the maximum elapsed time and greater than 0 (this allows for
        % the first call not to be counted twice)
        TimeOk = encounterSub.ElapsedTime<max_gap &...
            encounterSub.ElapsedTime>0;
        
        % Create bool for the similarity scores above the simulation
        % threshold and get their index within the sub encounter
        SimThreshOk = encounterSub.SimScore>=simthresh;
        goodIDs = find(TimeOk.*SimThreshOk);
        
        % If there are no calls in the encounter that meet the similarity
        % and elapsed time thresholds then start a new sub encounter or
        % move onto the next acoustic encounter
        if isempty(goodIDs)
            
            % Update the cluster id with the true (origional) indeces of
            % the call (calls were Clustered equals 1)
            chain(clusterID).index = ...
                encounterSub.TrueIDs(encounterSub.Clustered==1);
            
            % Remove the call(s) from the sub encounter from the acoustic
            % encounter
            encounterSub(encounterSub.Clustered==1,:)=[];
            
            % If, after removing the clustered calls from the acoustic
            % encounter, there are calls remaining move the iterater to the
            % first unclustered call and initialize the first call to the
            % new cluster
            if height(encounterSub)>0
                encounterSub.ElapsedTime = encounterSub.ArrivalTimes-...
                    encounterSub.ArrivalTimes(1);
                
                % Update the Simiarity Scores
                encounterSub.SimScore = simmat(...
                    encounterSub.TrueIDs(1),...
                    encounterSub.TrueIDs(1:end))';
                
                
                encounterSub.Clustered(1)=1;
            end
            idxInEnc =1;
            clusterID=clusterID+1;
        
        
        % If there were calls remaining in the acoustic encounter, move the
        % iterator to the next most similar call (above threshold) and add
        % it to the acoustic encounter
        else
            
            
            % Find which of the 'goodID' had the highest similarity score
            [~, maxIDx]=max(encounterSub.SimScore(goodIDs));
            
            % get the index in the subencounter of the good id and add it
            % to the cluster
            idxInEnc = goodIDs(maxIDx);
            encounterSub.Clustered(idxInEnc) =1;
            
            % Recalculate the elapsed time from the new call in the cluster
            % to the remaining calls in the acoustic encounter
            encounterSub.ElapsedTime = encounterSub.ArrivalTimes-...
                encounterSub.ArrivalTimes(idxInEnc);
            
            % Get the similarity scores for the next call in the acoustic
            % encouter to the rest of the calls in the acoustic encounter
            simVals = [
                zeros(idxInEnc-1,1)
                simmat(encounterSub.TrueIDs(idxInEnc),...
                       encounterSub.TrueIDs(idxInEnc:end))'
                ];
            
            % Update the similarity scores to reflect the new call index
            % within the acoustic encounter
            encounterSub.SimScore = simVals;
            
        end
        
    end
    
    
    
end

 
    
end





