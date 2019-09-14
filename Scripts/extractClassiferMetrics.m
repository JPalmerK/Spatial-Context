function [prct_improvement, methCorrect, methUnaided]= extractClassiferMetrics(perf)    
% Extrac the classifier performmance metrics from the structure
% Error rate reduction 
    try
    mod = struct2table(perf.Perf);
    catch me
        mod = struct2table(perf);
    end
    
        
    
    nCorrect = length(find(mod.ClusterScore> 0 & mod.TrueSpp ==1)) +...
        length(find(mod.ClusterScore<=0 & mod.TrueSpp ==0));
    methCorrect =nCorrect/height(mod);
    
    methUnaided = (length(find(mod.Score>= .5 & mod.TrueSpp ==1))...
        + length(find(mod.Score<.5 & mod.TrueSpp ==0)))/height(mod);
    
    error_unaided = 1-methUnaided;
    error_cluster = 1-methCorrect;
    prct_improvement = (error_unaided- error_cluster)/error_cluster;
   
end