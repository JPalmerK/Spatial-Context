function at = updateArrTableGPL(obj)
% Update the arrival table using GPL detections

% Parent hydrophone
parent = obj.array_struct.master;
child_hyd = obj.array_struct.slave;

ParentArrival = obj.localize_struct.hyd(parent).rtimes'./obj.fs;
TDOA =obj.localize_struct.hyd(parent).delays;
% Arrival time in seconds for all hydrophones based on delays
ArrivalSec =zeros(size(TDOA,1), (length(obj.hydrophone_struct)))/0;
ArrivalSec(:,parent) =ParentArrival;
CrossScores =zeros(size(TDOA,1), (length(obj.hydrophone_struct)))/0;

for ii=1:length(child_hyd)
    ArrivalSec(:,child_hyd(ii)) = (ParentArrival+TDOA(:,ii)).*...
        (obj.localize_struct.hyd(parent).cross_score(:,ii)./...
        obj.localize_struct.hyd(parent).cross_score(:,ii));
    CrossScores(:,child_hyd(ii)) = obj.localize_struct.hyd(parent).cross_score(:,ii);
end


x=squeeze(obj.localize_struct.hyd(parent).coordinates(end,1,:));
x(:,2)=squeeze(obj.localize_struct.hyd(parent).coordinates(end,2,:));

% Index of the detection
%idx = obj.localize_struct.hyd(parent).dex'+1;

idx = obj.localize_struct.hyd(parent).dex';

%
% Create the arrivals table (at) from the localize structure
at = table(...
    ArrivalSec,...
    x,...
    CrossScores,...
    TDOA,...
    idx,...
    'VariableNames',{'ArrivalSec', 'Location','CrossScore', 'TDOA', 'dex'});
at.ID = zeros(height(at),1)/0;
% Clear out any detections greater than Xkm from the receiver




end
