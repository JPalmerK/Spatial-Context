%% Function for applying simulated detector/classifier
function perfStructure = estClassifierPerf(obj)
% Compare the classification performance of the system without
% clustering to the system where all calls in a cluster are
% classified together


truthAndpreds = obj.truthTable;
truthAndpreds.PredClust = obj.Cluster_id;

cluster_ids = gather(unique(truthAndpreds.PredClust));
% For each cluster create the classificcation balues
for ii=1:length(cluster_ids)
    
    % get the ichain index value
    idx =find(truthAndpreds.PredClust==cluster_ids(ii));
    
    % Gt the classifier scores
    scores = truthAndpreds.Score(idx);
    
    
    % Likelihood ratio for the cluster
    LR = log(prod(scores./(1-scores)));
    
    % Use Geometric mean instead
    truthAndpreds.GeoMean(idx) = geomean(scores);
    
    % Fill in the likelihood ratio
    truthAndpreds.ClusterScore(idx) = LR;
    
    %Simple Voting
    pos_votes = sum(scores>0.5);
    neg_votes=sum(scores<=.5);
    
    if pos_votes>neg_votes
        truthAndpreds.Voting(idx) = ones(size(scores));
    elseif pos_votes<neg_votes
        truthAndpreds.Voting(idx) = zeros(size(scores));
    elseif pos_votes==neg_votes
        truthAndpreds.Voting(idx) = zeros(size(scores))+binornd(1,.5);
    end
    
end

obj.truthAndpreds = truthAndpreds;
perfStructure = table2struct(truthAndpreds);


end
