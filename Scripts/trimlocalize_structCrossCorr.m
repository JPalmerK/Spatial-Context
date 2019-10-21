function localize_struct_temp = trimlocalize_structCrossCorr(localize_struct, hyd, thresh)

localize_struct_temp = localize_struct;


% trim the scores

for ii = 1: length(localize_struct_temp.hyd)
    
    if ~isempty(localize_struct_temp.hyd(ii).delays)
        
        % Figure out where perfect cross correlation
        
        lowCrossVals = (localize_struct_temp.hyd(ii).cross_score)<thresh |...
            isnan(localize_struct_temp.hyd(ii).cross_score);
        
        localize_struct_temp.hyd(ii).cross_score(lowCrossVals) = nan;

        localize_struct_temp.hyd(ii).delays(lowCrossVals) = nan;
        
        % index of rows with TDOA's only
        k2 = find(sum(isnan(localize_struct_temp.hyd(ii).cross_score),2) ...
            < length(localize_struct_temp.hyd));

        % Trim calls
        localize_struct_temp.hyd(ii).score = localize_struct_temp.hyd(ii).score(:,k2);
        
        % Trim times
        localize_struct_temp.hyd(ii).rtimes = localize_struct_temp.hyd(ii).rtimes(:,k2);
        
        % Trim corrdinates
        localize_struct_temp.hyd(ii).coordinates = localize_struct_temp.hyd(ii).coordinates(:,:,k2);
        
        % trim dex
        localize_struct_temp.hyd(ii).dex = localize_struct_temp.hyd(ii).dex(k2);
        
        % trim coord time
        localize_struct_temp.hyd(ii).coord_time = localize_struct_temp.hyd(ii).coord_time(k2,:);
        
        % trim cross correlation score
        localize_struct_temp.hyd(ii).cross_score = localize_struct_temp.hyd(ii).cross_score(k2,:);
        
        % and delays (not dealing with CC matrix atm)
        localize_struct_temp.hyd(ii).delays = localize_struct_temp.hyd(ii).delays(k2,:);
        
        if isfield(localize_struct_temp.hyd(ii), 'detectorScore')
            localize_struct_temp.hyd(ii).detectorScore=...
                localize_struct_temp.hyd(ii).detectorScore(k2);
        end
    end
    
end
end