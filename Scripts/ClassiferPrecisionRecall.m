function [Precision, Recall]= ClassiferPrecisionRecall(perf)    
% Return the precision and recall for the clustered calls and the 
    try
    mod = struct2table(perf.Perf);
    catch me
        mod = struct2table(perf);
    end
    
    
    
    
    
    % Precision and Recall of the clustering algorithim
    nSpp_1 = sum(mod.TrueSpp);
    nCorrectSpp_1 = length(find((mod.TrueSpp ==1 & mod.GeoMean> .85)));
   
    Precision_clust = nCorrectSpp_1/ height(mod);
    Recall_clust = nCorrectSpp_1/nSpp_1;
    
    % Precision and Recall of the baseline classifier
    nCorrectSpp_1_unaided = sum(mod.Score> 0.85 & mod.TrueSpp ==1);
   
    Precision_unaided = nCorrectSpp_1_unaided/ height(mod);
    Recall_unaided = nCorrectSpp_1_unaided/nSpp_1;
    
    Precision = [Precision_clust Precision_unaided];
    Recall = [Recall_clust Recall_unaided];
    
    
end