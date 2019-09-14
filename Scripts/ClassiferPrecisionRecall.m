function [Precision, Recall]= ClassiferPrecisionRecall(perf)    
% Return the precision and recall for the clustered calls and the 
    try
    mod = struct2table(perf.Perf);
    catch me
        mod = struct2table(perf);
    end
    
    
    scorethreshs = [0.2:.05:.95];
    Precision_unaided =zeros(size(scorethreshs));
    Recall_unaided=zeros(size(scorethreshs));
    Precision_clust=zeros(size(scorethreshs));
    Recall_clust=zeros(size(scorethreshs));
    
    for ii = 1:length(scorethreshs)
        scorethresh = scorethreshs(ii);
    
    
    
    % Precision and Recall of the clustering algorithim
    nSpp_1 = sum(mod.TrueSpp);
    nCorrectSpp_1 = length(find(( mod.GeoMean> scorethresh & mod.TrueSpp ==1)));
   
    Precision_clust(ii) = nCorrectSpp_1/ height(mod);
    Recall_clust(ii) = nCorrectSpp_1/nSpp_1;
    
    % Precision and Recall of the baseline classifier
    nCorrectSpp_1_unaided = sum(mod.Score> scorethresh & mod.TrueSpp ==1);
   
    Precision_unaided(ii) = nCorrectSpp_1_unaided/ height(mod);
    Recall_unaided(ii) = nCorrectSpp_1_unaided/nSpp_1;

    end
    
    
end