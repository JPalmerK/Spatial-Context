function chain = updateChainsEncounterFirst(obj)
% Use the similarity matrix to update the cluster chains


% Use the similarity matrix to update the cluster chains


% Ladder Linkages rebuild
simmat= gather(obj.Sim_mat);
arraivalArray = gather(obj.arrivalArray);
max_gap = obj.maxEltTime;
simthresh = obj.cutoff;
ids = 1:size(arraivalArray,1);
timeDiffs = [0; diff(arraivalArray(:,1))];
acousticEncounterBreaks = [1; (find(timeDiffs>max_gap))];


dataTable = table();
dataTable.timeDiffs = timeDiffs;
dataTable.TrueIDs = ids';
dataTable.AcousticEncounters = ids'./ids';
dataTable.Clustered = ids'-ids';
dataTable.ArrivalTimes = arraivalArray(:,1);
dataTable.TrueCluster = arraivalArray(:,end);


for ii=1:length(acousticEncounterBreaks)
    
    dataTable.AcousticEncounters(acousticEncounterBreaks(ii):end) =ii;
    
end

clusterID =1;
chain =struct();
chainval = [];

% current index into the encouter
idx =1;


for ii=1:length(unique(dataTable.AcousticEncounters))
    
    % pull out the first acoustic encounter
    encounterSub = dataTable(dataTable.AcousticEncounters==ii,:);
    encounterSub.SimScore=encounterSub.Clustered;
    encounterSub.ElapsedTime = encounterSub.ArrivalTimes- ...
        encounterSub.ArrivalTimes(idx);
    
    % Get the similartiy scores for first detection in the dataset
    encounterSub.SimScore(idx:end) = simmat(encounterSub.TrueIDs(idx),...
        encounterSub.TrueIDs(idx):...
        encounterSub.TrueIDs(end))';
    
    encounterSub.Clustered(1)=1;
    
    % Index into the acoustic encounter
    idx =1;
    
    while height(encounterSub)>0
        % Identify all chains in the cluster
        TimeOk = encounterSub.ElapsedTime<max_gap & encounterSub.ElapsedTime>0;
        SimThreshOk = encounterSub.SimScore>=simthresh;
        goodIDs = find(TimeOk.*SimThreshOk);
        
        % if none are available end the cluster
        if isempty(goodIDs)
            chain(clusterID).index = encounterSub.TrueIDs(encounterSub.Clustered==1);
            encounterSub(encounterSub.Clustered==1,:)=[];
            
            % Update the similarity scores and delta times
            % update the elapsed time
            
            if height(encounterSub)>0
                encounterSub.ElapsedTime = encounterSub.ArrivalTimes-...
                    encounterSub.ArrivalTimes(1);
                
                % Update the Simiarity Scores
                encounterSub.SimScore = simmat(...
                    encounterSub.TrueIDs(1),...
                    encounterSub.TrueIDs(1:end))';
                
                
                encounterSub.Clustered(1)=1;
            end
            idx =1;
            clusterID=clusterID+1;
        else
            % find the where the maximum value of the comparison is and move
            % there
            
            [~, maxIDx]=max(encounterSub.SimScore(goodIDs));
            
            idx = goodIDs(maxIDx);
            encounterSub.Clustered(idx) =1;
            
            % Update the similarity scores and delta times
            % update the elapsed time
            encounterSub.ElapsedTime = encounterSub.ArrivalTimes-...
                encounterSub.ArrivalTimes(idx);
            
            % Update the Simiarity Scores
            simVals = [...
                zeros(idx-1,1);
                simmat( encounterSub.TrueIDs(idx),encounterSub.TrueIDs(idx:end))'];
            
            encounterSub.SimScore = simVals;
            
        end
        
    end
    
    
    
end

 
    
end





