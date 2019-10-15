function [GPL_struct] = GPL_kernel(GPL_struct, start, parm)

for k=1:length(start)
    GPL_struct(k).calltype=0; %Initially set all values to zero
end

[calltype_nums,calltype_index,~]=unique([parm.kernels.calltype]); % The unique call types


for calltype = 1:length(calltype_nums) % Loop through Call Types
    
    threshold=parm.kernels(calltype_index(calltype)).threshold; % One threshold per call type
    calltype_locs=find([parm.kernels.calltype]==calltype); % Indeces for this call type
    
    maxC=zeros(length(start),length(calltype_locs)); %allocate space to save cross correlations
    
    for kernel_num = 1:length(calltype_locs) % Loop through Kernels for that call type
        
        kernel=parm.kernels(calltype_locs(kernel_num)).call;
        kernel=kernel.^parm.kernels(kernel_num).power; % Raise kernel to specified power
        kernel_squared=kernel.^2; % Square kernel values
        kernel_sum_squares=sum(sum(kernel_squared)); % Sum all values
        sqrt_kernel_sum_squares=sqrt(kernel_sum_squares); %Take the square root of all kernel values
        
        norm_kernel=kernel./sqrt_kernel_sum_squares; % Normalize Kernel
        
        
        for k = 1:length(start) % Loop through Detections
            
            cm=GPL_full('cm',k,GPL_struct); % get cm matrix for detection
            cm=cm.^parm.kernels(kernel_num).power; %Raise cm to power
            cm_squared=cm.^2;
            cm_sum_squares=sum(sum(cm_squared));
            sqrt_cm_sum_squares=sqrt(cm_sum_squares);
            
            norm_cm=cm./sqrt_cm_sum_squares; %Normalize cm matrix
            
            zeropad=zeros(size(norm_kernel));
            
            norm_cm_padded=[zeropad,norm_cm,zeropad]; %NOTE: the number of freq bins must be the same in both kernel and cm
            
            % Cross correlate cm with kernel
            startcol=1;
            endcol=startcol+size(norm_kernel,2)-1;
            
            num_crosscorr=size(norm_cm_padded,2)-size(norm_kernel,2)+1; %total number of cross correlations to do
            C=zeros(1,num_crosscorr); %Allocate space for cross-corr values
            
            for i=1:num_crosscorr
                C(i)=sum(sum(norm_kernel.*norm_cm_padded(:,startcol:endcol)));%Cross-correlation
                startcol=startcol+1;
                endcol=endcol+1;
            end
            
            maxC(k,kernel_num)=max(C); %max cross correlation value
            
        end % looping through detections
        
    end %looping through kernels for call type
    
    maxC_percm=max(maxC,[],2); %Gives the maxium cross-corr value for each cm
    
    for k=1:length(start)
    GPL_struct(k).(strcat('calltype_',num2str(calltype_nums(calltype)),'_score'))=maxC_percm(k);%Save the max cross-corr value in output
    end
    
    truecall=find(maxC_percm>=threshold);%Indices of detections that exceed threshold
    
    for j=1:length(truecall)
        if GPL_struct(truecall(j)).calltype==0
            GPL_struct(truecall(j)).calltype=calltype_nums(calltype); %Add values for call type into output structure
        else
            GPL_struct(truecall(j)).calltype=max(calltype_nums)+1; %if this call was already marked as another call type, mark as next available number
        end
    end
    
end %looping through call types



end