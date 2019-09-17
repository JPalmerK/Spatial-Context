function localize_struct_temp = trimlocalize_struct(localize_struct, hyd, corrThresh, parentId)

localize_struct_temp = localize_struct;

array_id=parentId; % center hydrophone

% trim the scores

% Get the id's of the calls that were detected by two or mor hydrophone and
% with cross correlation scores above 

scores = cat(1,hyd(parentId).detection.calls.calltype_1_score);
UpcallIDX = find(scores>corrThresh);

% Now get the call id's represented by the calls for good matches. Localize
% struct contains the calls that were detected on two or more channels
UpcallCallIDs = intersect(UpcallIDX, localize_struct_temp.hyd(array_id).dex);


k2 = find(ismember(localize_struct_temp.hyd(array_id).dex,UpcallCallIDs));


% Trim calls
localize_struct_temp.hyd(array_id).score = localize_struct_temp.hyd(array_id).score(:,k2);

% Trim times
localize_struct_temp.hyd(array_id).rtimes = localize_struct_temp.hyd(array_id).rtimes(:,k2);

% Trim corrdinates
localize_struct_temp.hyd(array_id).coordinates = localize_struct_temp.hyd(array_id).coordinates(:,:,k2);

% trim dex
localize_struct_temp.hyd(array_id).dex = localize_struct_temp.hyd(array_id).dex(k2);

% trim coord time
localize_struct_temp.hyd(array_id).coord_time = localize_struct_temp.hyd(array_id).coord_time(k2,:);

% trim cross correlation score
localize_struct_temp.hyd(array_id).cross_score = localize_struct_temp.hyd(array_id).cross_score(k2,:);

% and delays (not dealing with CC matrix atm)
localize_struct_temp.hyd(array_id).delays = localize_struct_temp.hyd(array_id).delays(k2,:);
end