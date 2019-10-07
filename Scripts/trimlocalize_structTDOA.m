function localize_struct_temp = trimlocalize_structTDOA(localize_struct, hyd, driftSec)

localize_struct_temp = localize_struct;


% trim the scores

for ii = 1: length(localize_struct_temp.hyd)
    
    if ~isempty(localize_struct.hyd(ii).delays)
        

        for jj = 2:length(localize_struct_temp.hyd(ii).array.toa_diff)
            
            % Minimum allowable TDOA
            minTDOA = min(min(...
                localize_struct_temp.hyd(ii).array.toa_diff{jj})); 
           
            
            % calls under the tdoa
            badIdx = localize_struct_temp.hyd(ii).delays(:,jj-1)<=...
                (minTDOA -   driftSec); 
            localize_struct_temp.hyd(ii).delays(badIdx,jj-1) = nan;
            
            maxTDOA = max(max(...
                localize_struct_temp.hyd(ii).array.toa_diff{jj}));
            
            badIdx = localize_struct_temp.hyd(ii).delays(:,jj-1)>=...
                (maxTDOA +  driftSec); 
            localize_struct_temp.hyd(ii).delays(badIdx,jj-1) = nan;
            
        end
        
        % index of rows with TDOA's only
        k2 = find(sum(~isnan(localize_struct_temp.hyd(ii).delays),2)>0);

        % Trim calls
        localize_struct_temp.hyd(ii).score = ...
            localize_struct_temp.hyd(ii).score(:,k2);
        
        % Trim times
        localize_struct_temp.hyd(ii).rtimes =...
            localize_struct_temp.hyd(ii).rtimes(:,k2);
        
        % Trim corrdinates
        localize_struct_temp.hyd(ii).coordinates =...
            localize_struct_temp.hyd(ii).coordinates(:,:,k2);
        
        % trim dex
        localize_struct_temp.hyd(ii).dex = ...
            localize_struct_temp.hyd(ii).dex(k2);
        
        % trim coord time
        localize_struct_temp.hyd(ii).coord_time = ...
            localize_struct_temp.hyd(ii).coord_time(k2,:);
        
        % trim cross correlation score
        localize_struct_temp.hyd(ii).cross_score = ...
            localize_struct_temp.hyd(ii).cross_score(k2,:);
        
        % and delays (not dealing with CC matrix atm)
        localize_struct_temp.hyd(ii).delays = ...
            localize_struct_temp.hyd(ii).delays(k2,:);
    end
    
end
end