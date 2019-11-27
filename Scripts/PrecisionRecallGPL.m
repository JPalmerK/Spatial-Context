function [Clustered, Detector, RavenTable, NMI] = ...
    PrecisionRecallGPL(examp, parent, hyd)

% Input -
% Examp: simulation class object pre-processed for similarity matrix
% Truth: Validated Raven selection table
% Threshold: correlation threshold



localize_struct = examp.localize_struct;
arivalArr = examp.arrivalArray;
cluster_ids = unique(examp.Cluster_id);

%% Create a raven table from the detections

% Convert arrival time to seconds
start_times =localize_struct.hyd(parent).rtimes'/examp.fs;



% Lat and Lon
lat = squeeze(localize_struct.hyd(parent).coordinates(end,1,:));
lon = squeeze(localize_struct.hyd(parent).coordinates(end,2,:));

scores = localize_struct.hyd(parent).detectorScore;
scores = reshape(scores, [length(scores),1]);
ClusterId = examp.Cluster_id;

% Validation info 1 is target, others are not. In graywhale dataset call
% type 2 is fin and call type 0 is noise. 
Pruned = localize_struct.hyd(parent).pruned;
Pruned = reshape(Pruned, [length(Pruned),1]);
% Index of the call
dex = localize_struct.hyd(parent).dex';
RavenTable = table(...
    start_times,...
    lat,...
    lon,...
   scores,...
    ClusterId,...
    Pruned,...
    dex,...
    'VariableNames',{ ...
    'BeginS', 'Lat', 'Lon','Scores','ClusterId', 'pruned','dex'});


% Avoid division by 0
RavenTable.Scores= RavenTable.Scores+0;
RavenTable.voting= zeros([height(RavenTable),1])/0;
ClusterEntropy = zeros(size(cluster_ids));


for jj = 1:length(cluster_ids)
    
    clus = cluster_ids(jj);
    
    corrScores = RavenTable.Scores(RavenTable.ClusterId==clus);
    corrScores = corrScores(~isnan(corrScores));
    
    LR =  log(prod(corrScores./(1.00001-corrScores)))/length(corrScores);

    RavenTable.voting(RavenTable.ClusterId==clus)=LR;
    
    % Entropy of cluster labels
    Hc(jj) = sum(RavenTable.ClusterId==clus)/height(RavenTable);
    
    % Calculate the entropy of class labels within each cluster
    prunedvals= ( RavenTable.pruned(RavenTable.ClusterId==clus));
    prunedIds = unique(prunedvals);
    
    % Determine the proprtion of the cluster made up by each class within
    % the cluster
    for kk=1:length(prunedIds)
        callprops(kk)= sum(prunedvals==prunedIds(kk))/length(prunedvals);
    end
    
    
    % Conditional Class Entropy for each cluster
    Hyc(jj) =  -Hc(jj) * sum(callprops.*log2(callprops));
    
    

end

RavenTable.voting(RavenTable.voting==Inf) = ...
    max(RavenTable.voting(isfinite(RavenTable.voting)));

% Entropy of class labels
Y = unique(RavenTable.pruned);
Hy =zeros(size(Y));

% Class Entropy- get the proportion of labels of each class in the dataset
for ii=1:length(Y)
    Hy(ii)= sum(RavenTable.pruned==Y(ii))/height(RavenTable);  
end

% Entropy of class labels
Hy= sum(-Hy.*log2(Hy));


% Entropy of class labels
Hc =sum(-Hc.*log2(Hc));




% Mutual information Hy -H(Y|C) - conditional entropy for each class
IYC = Hy-sum(Hyc);

% Normalized  Mutual Information 
NMI = ((2*IYC)/(Hy +Hc));




% threhsold scores for calculating precision recall
RavenTable.LRscores =  log(RavenTable.Scores./(1.00001-RavenTable.Scores));
    

scores = unique(RavenTable.LRscores);
scores = unique(round(scores,2));






%Pre allocate output
RecallDetOnly = zeros(1, length(scores)-1);
PrecisionDetOnly = RecallDetOnly;
RecallClustered= RecallDetOnly;
PrecisionClustered = RecallDetOnly;
ErrorRateClustered = RecallDetOnly;
ErrorRateDetector = RecallDetOnly;

for jj = 1:length(scores)
    
    % Detector only Precision Recall
    TP(jj) = sum((RavenTable.LRscores >= scores(jj))&...
        (RavenTable.pruned==1));
    
    
    % Labeled right whale but were not right whale
    FP(jj) = sum((RavenTable.LRscores >= scores(jj)).*...
        (RavenTable.pruned ~=1));
    
    % Not labeled right whale but were right whale
    FN(jj) = sum((RavenTable.LRscores < scores(jj)).*...
        (RavenTable.pruned==1));
    
    % Number of detections that were labeled RW and were in agreement with the
    % truth value
    TPv(jj) = sum((RavenTable.voting >= (scores(jj))).*...
        (RavenTable.pruned==1));
    
    % Labeled right whale but were not right whale
    FPv(jj) = sum((RavenTable.voting >= (scores(jj))).*...
        (RavenTable.pruned~=1));
    
    % Not labeled right whale but were right whale
    FNv(jj) = sum((RavenTable.voting <(scores(jj))).*....
        (RavenTable.pruned==1));
    
end


RecallDetOnly = TP./(TP+FN);
PrecisionDetOnly = TP./(TP+FP);



% Add detections that were trimmed via thresholding to the fn's
RecallClustered= TPv./(TPv+FNv);
PrecisionClustered = TPv./(TPv+FPv);


Clustered= struct();
Clustered.Precision = PrecisionClustered;
Clustered.Recall =RecallClustered;


Detector= struct();
Detector.Precision = PrecisionDetOnly;
Detector.Recall =RecallDetOnly;





end












