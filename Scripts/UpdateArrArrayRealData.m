function array = UpdateArrArray(simStruct)
% This funcion creates the array of arrivals, if the
% arrivalsMiss parameter is set to 1, then where multiple calls
% fall within the expected TDOA range, the call with the
% arrival TDOA pair nearest to TDOA of 0 is chosen.
% This selection is somewhat arbritrary as per design.


% Array containing the parent and children hydrophones

% Get ID vlue, if simulated data then ID if GPL then dex
% reference
if isfield(simStruct, {'calls'}) && ~isempty(simStruct.calls)
    disp('GPL Calls Detected, array ID variable being filled by dex')
    IDval = simStruct.arrivalTable.dex+1; % Bug in GPL code leave this in until recieved fix from Tyler XXX
else
    IDval = simStruct.arrivalTable.ID;
end


array = [simStruct.arrivalTable.ArrivalSec(:,...
    [simStruct.array_struct.master,...
    simStruct.array_struct.slave(simStruct.child_idx)])...
    simStruct.arrivalTable.Location ...
    IDval];


% Remove calls that weren't detected on the parent hydrophone
array = array(~isnan(array(:,1)),:);

% At least 1 tdoa value needed, therefore remove any calls
% (master) where there aren't at least two arrivals
array = array(sum(~isnan(array(:,1:end)),2)>= 2,:);




% Sort by start time on the parent hydrophone
[~,sortidx] = sort(array(:,1));
array = (array(sortidx,:));





end
