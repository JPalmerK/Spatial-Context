function [prct_improvement, methCorrect, methUnaided]= extractClassiferMetrics(perf)    
% Extrac the classifier performmance metrics from the structure

    mod = struct2table(perf.perf);
    nCorrect = length(find(mod.ClusterScore> 0 & mod.TrueSpp ==1)) +...
        length(find(mod.ClusterScore<=0 & mod.TrueSpp ==0));
    methCorrect =nCorrect/height(mod);
    methUnaided = (length(find(mod.Score>= .5 & mod.TrueSpp ==1))...
        + length(find(mod.Score<.5 & mod.TrueSpp ==0)))/height(mod);
    prct_improvement = -(methUnaided- methCorrect)/methUnaided;
   
end