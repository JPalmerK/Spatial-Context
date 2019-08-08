function localize_struct = trimlocalize_struct(localize_struct, lsqScore, parentId)

array_id=parentId; % center hydrophone

% trim the scores
score=localize_struct.hyd(array_id).score(5,:);
[~, k2]= find(score < lsqScore);

% Trim calls
localize_struct.hyd(array_id).score = localize_struct.hyd(array_id).score(:,k2);

% Trim times
localize_struct.hyd(array_id).rtimes = localize_struct.hyd(array_id).rtimes(:,k2);

% Trim corrdinates
localize_struct.hyd(array_id).coordinates = localize_struct.hyd(array_id).coordinates(:,:,k2);

% trim dex
localize_struct.hyd(array_id).dex = localize_struct.hyd(array_id).dex(k2);

% trim coord time
localize_struct.hyd(array_id).coord_time = localize_struct.hyd(array_id).coord_time(k2,:);

% trim cross correlation score
localize_struct.hyd(array_id).cross_score = localize_struct.hyd(array_id).cross_score(k2,:);

% and delays (not dealing with CC matrix atm)
localize_struct.hyd(array_id).delays = localize_struct.hyd(array_id).delays(k2,:);
end