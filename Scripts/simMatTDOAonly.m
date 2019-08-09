function Sim_mat = simMatTDOAonly(simStruct)

% Create the simulation matrix using TDOA values only

if isempty(simStruct.TDOA_vals)
    disp(['Updating TDOA values'])
    UpdateTDOA(simStruct);
end

try
    
    Sim_mat = zeros(length(simStruct.arrivalArray(:,end)), 'gpuArray')/0;
catch
    Sim_mat = zeros(length(simStruct.arrivalArray(:,end)))/0;
end

for ii =1:(length(simStruct.arrivalArray)-1)
    
    
    tdoa_orig = simStruct.TDOA_vals(ii,:);
    
   
    
    % index of all calls within the elapsed time
    nextTimes = simStruct.arrivalArray(ii:end,1);
    time_diffs =  nextTimes- simStruct.arrivalArray(ii,1);
    
    % Identify the calls within the the acoustic encounter
    acousticEncIdx = find(diff(nextTimes)>= simStruct.maxEltTime,1);
    
    time_diffs = time_diffs(1:acousticEncIdx-1);
    
    % Get the TDOA values
    TDOA_next = simStruct.TDOA_vals((ii+(acousticEncIdx-1)),:);
    
    
    
    % For each hydrophone pair calculate likelihood values
    deltaTDOALklhd =[];
    
    
    for jj=1:size(TDOA_next,2)
        
        mu = zeros(length(time_diffs),1);
        sigmaSwim = 2* sqrt((time_diffs * (simStruct.s)/simStruct.c).^2);
        x =(tdoa_orig(jj) - TDOA_next(:,jj)); %values
        
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
    Sim_mat(ii,ii) = 1;
    Sim_mat(ii, ii:ii+length(simValues)-1) = simValues;
    Sim_mat(ii:ii+length(simValues)-1,ii) = simValues;
end



end