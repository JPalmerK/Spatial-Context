function chain = updateChainsEncounterFirst(obj)
% Use the similarity matrix to update the cluster chains


% Ladder Linkages rebuild
simmat= gather(obj.Sim_mat);
arraivalArray = gather(obj.arrivalArray);
max_gap = obj.maxEltTime;
simthresh = obj.cutoff;
ids = 1:length(arraivalArray);
timeDiffs = [0; diff(arraivalArray(:,1))];
acousticEncounterBreaks = unique([1; (find(timeDiffs>max_gap))-1; length(ids)]);


dataTable = table();
dataTable.timeDiffs = timeDiffs;
dataTable.TrueIDs = ids';
dataTable.AcousticEncounters = ids'./ids';
dataTable.Clustered = ids'-ids';
dataTable.ArrivalTimes = arraivalArray(:,1);
dataTable.TrueCluster = arraivalArray(:,end);


for ii=1:length(acousticEncounterBreaks)-1
    
    dataTable.AcousticEncounters(acousticEncounterBreaks(ii):...
        acousticEncounterBreaks(ii+1) ) =ii;
    
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
    idx =1;
    
    chainval=[];
    if height(encounterSub)> 1
        
        
        while height(encounterSub)>1
            
            
            
            % get the id of the call with the highest similarity score where the
            % maximum time gap hasn't been reached
            TimeOk = encounterSub.ElapsedTime<max_gap & encounterSub.ElapsedTime>0;
            SimThreshOk = encounterSub.SimScore>=simthresh;
            goodIDs = find(TimeOk.*SimThreshOk);
            
            %if good ID's is empty start a new cluster and move onto the next index
            if isempty(goodIDs)
                
                
                chainval = unique([chainval encounterSub.TrueIDs(idx)]);
                chain(clusterID).index =  chainval;
                chain(clusterID).n = length(chainval);
                encounterSub.Clustered(idx)=1;
                
                idx = idx +1;
                
                % if the index is beyond the encouter then we need to start over
                if idx>=height(encounterSub)
                    idx =2;
                end
                
                % update the elapsed time
                encounterSub.ElapsedTime = encounterSub.ArrivalTimes-...
                    encounterSub.ArrivalTimes(idx);
                
                % Update the Simiarity Scores
                encounterSub.SimScore(idx:end) = simmat(...
                    encounterSub.TrueIDs(idx),...
                    encounterSub.TrueIDs(idx:end))';
                
                
                
                % Remove the rows from the ones we just clustered
                encounterSub(encounterSub.Clustered==1,:) =[];
                clusterID=clusterID+1;
                
                encounterSub.Clustered(1)=1;
                chainval=[];
                
                idx =1;
                
                % otherwise calls within encounter that match similarity thresholds and
                % time thresholds
            elseif length(goodIDs)>1
                
                [~, maxIDx]=max(encounterSub.SimScore(goodIDs));
                chainval= [chainval,...
                    encounterSub.TrueIDs(goodIDs(maxIDx))];
                
                idx = goodIDs(maxIDx);
                
                % Update the Simiarity Scores
                encounterSub.SimScore(idx:end) = simmat(...
                    encounterSub.TrueIDs(idx),...
                    encounterSub.TrueIDs(idx:end))';
                
                encounterSub.ElapsedTime = encounterSub.ArrivalTimes-...
                    encounterSub.ArrivalTimes(idx);
                encounterSub.Clustered(idx)=1;
                
                
                
                
            elseif length(goodIDs)==1
                
                chainval = [chainval encounterSub.TrueIDs(goodIDs)];
                chain(clusterID).index =  chainval;
                chain(clusterID).n = length(chainval);
                encounterSub.Clustered(goodIDs)=1;
                idx = goodIDs;
                
                % if the index is beyond the encouter then we need to start over
                if idx>=height(encounterSub)
                    idx =2;
                else
                    idx = goodIDs+1;
                end
                
                encounterSub.ElapsedTime = encounterSub.ArrivalTimes-...
                    encounterSub.ArrivalTimes(idx);
                
                % Update the Simiarity Scores
                encounterSub.SimScore(idx:end) = simmat(...
                    encounterSub.TrueIDs(idx),...
                    encounterSub.TrueIDs(idx:end))';
                % Remove the rows from the ones we just clustered
                encounterSub(encounterSub.Clustered==1,:) =[];
                
                if height(encounterSub)>0
                encounterSub.Clustered(1)=1;
                end
                clusterID=clusterID+1;
                chainval=[];
                idx =1;
            end
            
        end
        
    else
        chain(clusterID).index =  encounterSub.TrueIDs(1);
        chain(clusterID).n = 1;
        clusterID=clusterID+1;
    end
    
    
    
    
end





