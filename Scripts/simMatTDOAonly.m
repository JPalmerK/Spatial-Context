function Sim_mat = simMatTDOAonly(simStruct)

% Create the simulation matrix using TDOA values only

Sim_mat = zeros(size(simStruct.arrivalArray))/0;

arrivalArray= (simStruct.arrivalArray);
for ii =1:(size(simStruct.TDOA_vals,1))
    
    
    tdoa_orig =simStruct.TDOA_vals(ii,:);
    
    
    % Figure out the number of time gaps within the maximum
    % allowed correlation time (time_cut)
    time_gaps = arrivalArray(ii:end, 1)-...
        arrivalArray(ii, 1);
    
    diff_times = diff(arrivalArray(ii:end, 1));
    
    % Find first big gap
    idx_end = find(diff_times>= simStruct.maxEltTime,1)-1;
    
    if isempty(idx_end)
        idx_end = length(time_gaps);
    end
    
    time_gaps = time_gaps(1:idx_end);
    TDOA_next = simStruct.TDOA_vals(ii:ii+idx_end-1,:);
    
    
    deltaTDOA = bsxfun(@minus, tdoa_orig,TDOA_next);
    elapsedTime = time_gaps;
    
    deltaTDOALklhd=zeros(size(deltaTDOA));
    for jj=1:size(deltaTDOA,2)
        
        mu = zeros(size(TDOA_next,1),1);
        sigmaSwim = 2* sqrt(simStruct.drift^2+ (elapsedTime * (simStruct.s)/simStruct.c).^2);
        x =deltaTDOA(:,jj); %values
        
        likelihood = normpdf(x,mu,sigmaSwim);
        
        % Normalizing factor
        LikelihoodNormFac= normpdf(0,0,sigmaSwim);
        NormLikelihood = likelihood./LikelihoodNormFac;
        
        % Normalized likelihood
        deltaTDOALklhd(:,jj)= NormLikelihood; % sigma
        % Create normalizing factor
   
    end
    
    simValues = nanmin(deltaTDOALklhd,[],2);
    
    % Take the minimum value and fill in the similarity matrix
    
    Sim_mat(ii, ii:ii+length(simValues)-1) = simValues;
    Sim_mat(ii:ii+length(simValues)-1,ii) = simValues;
    Sim_mat(ii,ii) = 1;
end



end