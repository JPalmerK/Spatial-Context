function localize_struct_temp = trimlocalize_struct(localize_struct, hyd, corrThresh)
% Trim the localize structure based on detector scores
localize_struct_temp = localize_struct;


% trim the scores

for ii = 1: length(localize_struct_temp.hyd)
    
    if ~isempty(localize_struct.hyd(ii).delays)
        
        % Get the id's of the calls that were detected by two or mor hydrophone and
        % with cross correlation scores above
        
        k2 = find(localize_struct_temp.hyd(ii).detectorScore>=corrThresh);
        
        
        
        % More comments needed!
        
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
        
        % and delays (not dealing with CC matrix atm)
        localize_struct_temp.hyd(ii).detectorScore = ...
            localize_struct_temp.hyd(ii).detectorScore(k2);
        
        
        if isfield(localize_struct_temp.hyd(ii), 'pruned')
            localize_struct_temp.hyd(ii).pruned=...
                localize_struct_temp.hyd(ii).pruned(k2);
        end
        
    end
    
end
end