% Function for returning the sensitiviery loop

function [ExpScoresMeth, UnAidedPerf] = runClassPerfLoop(simStruct,TimeThresh,SimThresh)
% Returns mxn matrix of adjusted rand scores where m is the length of the
% time threshold and n is the length of the similarity threshold
% Input - simulation object, TimeThresh, SimThres



% If only two inputs then baseline, assign simThresh 1 and skip update
% chains
if nargin ==2
    SimThresh =1;
    
    ExpScoresMeth = zeros(length(TimeThresh), length(SimThresh))/0;
    UnAidedPerf =  ExpScoresMeth;
    nAgents =ExpScoresMeth;
    
    
    simStruct.truthTable = createTruthTable(simStruct);
    simStruct.truthTable=createSpeciesPreds(simStruct);
    simStruct.Cluster_id = acEnc(simStruct);
    
    
    for jj = 1:length(SimThresh)
        simStruct.Cluster_id =[];
        simStruct.cutoff = SimThresh(jj);
        
        for ii = 1:length(TimeThresh)
            
            simStruct.Cluster_id =[];
            simStruct.maxEltTime=(TimeThresh(ii));
            
            simStruct.Cluster_id = acEnc(simStruct);
            
            simStruct.truthAndpreds=createSpeciesPreds(simStruct);
            perf = estClassifierPerf(simStruct);
            [deltaPerf unaided] = extractClassiferMetrics(perf);
            

            
            
            ExpScoresMeth(ii,1) = deltaPerf;
            UnAidedPerf(ii,1) =unaided;
            
        end
        
    end
    
    
else % Else there were three variables, so leave well enough alone

    idx =0;
    totit = length(SimThresh)*length(TimeThresh);
    
    
    ExpScoresMeth = zeros(length(TimeThresh), length(SimThresh))/0;
    UnAidedPerf =  ExpScoresMeth;
    % add the species predictions
    simStruct.truthTable=createSpeciesPreds(simStruct);
    
    for jj = 1:length(SimThresh)
        simStruct.Cluster_id =[];
        simStruct.cutoff = SimThresh(jj);
        
        for ii = 1:length(TimeThresh)
         
            simStruct.Cluster_id =[];
            simStruct.maxEltTime=(TimeThresh(ii));
            
            simStruct.chains =updateChainsEncounterFirst(simStruct);
            simStruct.Cluster_id= updateClusterID(simStruct);
           

            
            perf = estClassifierPerf(simStruct);
            [deltaPerf, methCorrect, methUnaided] = extractClassiferMetrics(perf);
            
            
%             [Precision, Recall]= ClassiferPrecisionRecall(perf);
%             
%             Precision_cluster(ii,jj) = Precision(1);
%             Recall_cluster(ii,jj) =Recall(1);
%             
%             Precision_unaided(ii,jj) =Precision(2);
%             Recall_unaided(ii,jj) = Recall(2);

            
            
            ExpScoresMeth(ii,jj) = methCorrect;
            UnAidedPerf(ii,jj) =methUnaided;
            
            
        end
        disp([num2str(idx) ' of ' num2str(totit) 'time/simthresh'])
    end
end

end




