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
            
            simStruct.Cluster_id = acEnc(simStruct);
            ExpScoresMeth(ii,1)  = getRand(simStruct);
            nAgents(ii,1) = length(unique(simStruct.Cluster_id));
            
        end
        
    end
    
    
else % Else there were three variables, so leave well enough alone
    try
        ExpScoresMeth = zeros(length(TimeThresh), length(SimThresh), 'gpuArray')/0;
    catch
        
        ExpScoresMeth = zeros(length(TimeThresh), length(SimThresh))/0;
    end
    nAgents =ExpScoresMeth;
    idx =0;
    totit = length(SimThresh)*length(TimeThresh);
    
    for jj = 1:length(SimThresh)
        simStruct_copy = simStruct;
        simStruct_copy.Cluster_id =[];
        simStruct_copy.cutoff = SimThresh(jj);
        
        for ii = 1:length(TimeThresh)
            simStructTmp = simStruct_copy;
            simStructTmp.Cluster_id =[];
            simStructTmp.maxEltTime=(TimeThresh(ii));
            
            simStructTmp.chains =updateChainsEncounterFirst(simStructTmp);
            simStructTmp.Cluster_id= updateClusterID(simStructTmp);
            
            ExpScoresMeth(ii,jj) =  getRand(simStructTmp);
            nAgents(ii,jj) = length(unique(simStructTmp.Cluster_id));
            
            
            idx=idx+1;
        end
        disp([num2str(idx) ' of ' num2str(totit) 'time/simthresh'])
    end
end





end




