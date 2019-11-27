function LocalizeOut = trimFalsePositivesGPLPrune(localize_struct, parent, hyd, truth, propRem)
% This function randomly removes some proportion of the false positive
% detections s defined by propRem



for ii=1: length(localize_struct.hyd)
    if ~isempty(localize_struct.hyd(ii).delays)
        
        truth_temp = truth_trimmed(truth_trimmed.Channel==parent,:);
        arrival_sec = localize_struct.hyd(ii).rtimes'/2000;
        arrivalMtlbDate = datenum('20090328', 'yyyymmdd')+arrival_sec/86400;
        sppMat=repmat({'unknown'}, [length(arrivalMtlbDate),1]);
        
        
        wiggleroom_days = .3/60/60/24;
        
        % step through and add the species where available
        for jj=1:length(arrivalMtlbDate)
            
            if height(truth_temp)>0
                % Logicals for finding calls overlapping with GPL detections
                gpl_time = arrivalMtlbDate(jj);
                
                linkedID = find(...
                    truth_temp.mend> (gpl_time-wiggleroom_days) &...
                    truth_temp.mstart< (gpl_time +wiggleroom_days));
                
                
                if length(linkedID)>=1
                    spp =truth_temp.Species(linkedID);
                    sppMat(jj) = spp(1);
                    truth_temp(linkedID,:)=[];
                end
            end
            
            
            
        end
        
        FP_ids = find(~strcmp(sppMat, 'rw'));
        
        % Number of calls to remove
        nDitch = ceil(length(FP_ids)*(propRem/100));
        
        ditchIdx = randsample(FP_ids, nDitch);
        
        
        localize_struct.hyd(ii).cross_score(ditchIdx,:)=[];
        localize_struct.hyd(ii).coord_time(ditchIdx,:)=[];
        localize_struct.hyd(ii).rtimes(ditchIdx)=[];
        localize_struct.hyd(ii).delays(ditchIdx,:)=[];
        localize_struct.hyd(ii).dex(ditchIdx)=[];
        localize_struct.hyd(ii).coordinates(:,:,ditchIdx)=[];
      
        localize_struct.hyd(ii).score(:,ditchIdx)=[];
        
        if isfield(localize_struct.hyd(ii), 'detectorScore')
         localize_struct.hyd(ii).detectorScore(ditchIdx)=[];
        end
        
        LocalizeOut=localize_struct;
        
    end
    
    
end