function Sim_mat = simMatTDOAonly(simStruct)

% Create the simulation matrix using TDOA values only

Sim_mat = zeros(size(simStruct.arrivalArray));

for ii =1:(length(simStruct.arrivalArray))
    
    % First TDOA
    tdoa_orig = simStruct.TDOA_vals(ii,:);
    
    
    % index of all calls within the elapsed time
    nextTimes = simStruct.arrivalArray(ii:end,1);
    elapsedTime =  nextTimes- nextTimes(1);
    
    % Identify the calls within the the acoustic encounter
    acousticEncIdx = find(diff(nextTimes)>= simStruct.maxEltTime,1)-1;
    elapsedTime = elapsedTime(1:acousticEncIdx);
    
    % Get the TDOA values
    nextTDOAidxStart = ii+1;
    nextTDOAidxEnd = ii+(acousticEncIdx)-1;
    TDOA_next = simStruct.TDOA_vals(nextTDOAidxStart:nextTDOAidxEnd,:);
    
    
    
    % For each hydrophone pair calculate likelihood values
    deltaTDOALklhd =[];
    
    
    for jj=1:size(TDOA_next,2)
        
        mu = zeros(length(elapsedTime),1);
        sigmaSwim = 2* sqrt((elapsedTime * (simStruct.s)/simStruct.c).^2);
        x =[0; (tdoa_orig(jj) - TDOA_next(:,jj))]; %values
        
        likelihood = normpdf(x,mu,sigmaSwim);
        
        % Normalizing factor
        LikelihoodNormFac= normpdf(0,0,sigmaSwim);
        NormLikelihood = likelihood./LikelihoodNormFac;
        
        % Normalized likelihood
        deltaTDOALklhd = [deltaTDOALklhd, NormLikelihood]; % sigma
        % Create normalizing factor
   
    end
    
    simValues = nanmin(deltaTDOALklhd,[],2);
    
    % Take the minimum value and fill in the similarity matrix
    
    Sim_mat(ii, ii:ii+length(simValues)-1) = simValues;
    Sim_mat(ii:ii+length(simValues)-1,ii) = simValues;
    Sim_mat(ii,ii) = 1;
end



end