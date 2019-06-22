% Function for returning the sensitiviery loop

function ExpScoresMeth = runSensitivtyLp(examp,TimeThresh,SimThresh)
% Returns mxn matrix of adjusted rand scores where m is the length of the
% time threshold and n is the length of the similarity threshold
% Input - simulation object, TimeThresh, SimThres


% If only two inputs then baseline, assign simThresh 1 and skip update
% chains
if nargin ==2
    SimThresh =1;
    
    ExpScoresMeth = zeros(length(TimeThresh), length(SimThresh))/0;
   
    
    for jj = 1:length(SimThresh)
        examp.Cluster_id =[];
        examp.cutoff = SimThresh(jj);
        
        for ii = 1:length(TimeThresh)
            
            examp.Cluster_id =[];
            examp.time_cut=(TimeThresh(ii));
            examp.toaOnlyCluster();
            examp.getRand();
            ExpScoresMeth(ii,1) = examp.AdjRand;
            
        end
        
    end
    
    
else % Else there were three variables, so leave well enough alone
    ExpScoresMeth = zeros(length(TimeThresh), length(SimThresh))/0;

    for jj = 1:length(SimThresh)
        examp.Cluster_id =[];
        examp.cutoff = SimThresh(jj);
        
        for ii = 1:length(TimeThresh)
            
            examp.Cluster_id =[];
            examp.time_cut=(TimeThresh(ii));
            examp.updateChains;
            examp.getRand();
            ExpScoresMeth(ii,jj) = examp.AdjRand;
            
        end
        
    end
end





end




