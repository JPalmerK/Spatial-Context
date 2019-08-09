% Function for returning the sensitiviery loop

function [ExpScoresMeth, nAgents] = runSensitivtyLp(simStruct,TimeThresh,SimThresh)
% Returns mxn matrix of adjusted rand scores where m is the length of the
% time threshold and n is the length of the similarity threshold
% Input - simulation object, TimeThresh, SimThres


% If only two inputs then baseline, assign simThresh 1 and skip update
% chains
if nargin ==2
    SimThresh =1;
    
    ExpScoresMeth = zeros(length(TimeThresh), length(SimThresh))/0;
    nAgents =ExpScoresMeth;
    
    for jj = 1:length(SimThresh)
        simStruct.Cluster_id =[];
        simStruct.cutoff = SimThresh(jj);
        
        for ii = 1:length(TimeThresh)
            
            simStruct.Cluster_id =[];
            simStruct.maxEltTime=(TimeThresh(ii));
            simStruct.toaOnlyCluster();
            simStruct.getRand();
            ExpScoresMeth(ii,1) = simStruct.AdjRand;
            nAgents(ii,1) = length(unique(simStruct.Cluster_id));
            
        end
        
    end
    
    
else % Else there were three variables, so leave well enough alone
    ExpScoresMeth = zeros(length(TimeThresh), length(SimThresh))/0;
    nAgents =ExpScoresMeth;
    idx =0;
    totit = length(SimThresh)*length(TimeThresh);
    for jj = 1:length(SimThresh)
        simStruct.Cluster_id =[];
        simStruct.cutoff = SimThresh(jj);
        
        for ii = 1:length(TimeThresh)
            
            simStruct.Cluster_id =[];
            simStruct.maxEltTime=(TimeThresh(ii));
            
            simStruct.chains =updateChainsEncounterFirst(simStruct);

            
            simStruct.Cluster_id= updateClusterID(simStruct);
            
            ExpScoresMeth(ii,jj) =  getRand(simStruct);
            nAgents(ii,jj) = length(unique(simStruct.Cluster_id));
            
            %disp([num2str(idx) ' of ' num2str(totit) 'time/simthresh'])
            idx=idx+1;
        end
        
    end
end





end




