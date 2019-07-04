% Create clusters of data from GPL localization struct
classdef clusterClassGPL <simulationClass
    properties
        CrossScore= .8
        localize_struct % Localize struct if using GPL
        ravenTable % Raven sound selection table
        limitTime % to limit the time over which to analyse the data fill here (seconds)
    end
    methods
        function UpdateArrTable(obj)
            
            if isempty(obj.localize_struct) & isempty(obj.ravenTable)
                disp('Raven Selection Table or GPL localization Structure needed')
                return
            elseif  isempty(obj.localize_struct) & ~isempty(obj.ravenTable)
                updateArrTableRaven(obj) % Not implemented
                disp('Ravent table access not yet implemented')
            elseif ~isempty(obj.localize_struct) & isempty(obj.ravenTable)
                updateArrTableGPL(obj)
            else ~isempty(obj.localize_struct) & ~isempty(obj.ravenTable)
                disp('Raven selection table OR GPL localize struct needed')
                return
            end
        end
        
        function updateArrTableGPL(obj)
            
            % Parent hydrophone
            parent = obj.array_struct.master;
            
            child_hyd = obj.array_struct.slave;
            
            ParentArrival = obj.localize_struct.hyd(parent).rtimes'./obj.fs;
            TDOA =obj.localize_struct.hyd(parent).delays;
            
            ArrivalSec =zeros(length(TDOA), (length(obj.hydrophone_struct)))/0;
            ArrivalSec(:,parent) =ParentArrival;
            CrossScores =zeros(length(TDOA), (length(obj.hydrophone_struct)))/0;
            
            for ii=1:length(child_hyd)
                ArrivalSec(:,child_hyd(ii)) = (ParentArrival+TDOA(:,ii)).*...
                    (obj.localize_struct.hyd(parent).cross_score(:,ii)./...
                    obj.localize_struct.hyd(parent).cross_score(:,ii));
                CrossScores(:,child_hyd(ii)) = obj.localize_struct.hyd(parent).cross_score(:,ii);
            end
            
            
            x=squeeze(obj.localize_struct.hyd(parent).coordinates(2,1,:));
            x(:,2)=squeeze(obj.localize_struct.hyd(parent).coordinates(2,2,:));
            %
            % Create the arrivals table (at) from the localize structure
            at = table(...
                ArrivalSec,...
                x,...
                CrossScores,...
                TDOA,...
                'VariableNames',{'ArrivalSec', 'Location','CrossScore', 'TDOA'});
            at.ID = zeros(height(at),1)/0;
            % Clear out any detections greater than Xkm from the receiver
            crossScore_bool = ones(size(at.CrossScore));
            crossScore_bool(at.CrossScore>obj.CrossScore) = nan;
            
            % Add an extra column for the parent hydropone
            at.CrossScore = at.CrossScore.*crossScore_bool;
            at.ArrivalSec = at.ArrivalSec.*crossScore_bool;
            
            % Remove rows wehre the cross correlation score isn't as good
            % as the minimum
            cross_scores = sum(crossScore_bool(:,child_hyd(obj.child_idx)),2);
            good_idx = logical(~isnan(cross_scores));
            
            at = at(good_idx,:);
            
            
            if ~isempty(obj.limitTime)
                good_idx = logical(cumsum(diff(ArrivalSec(:, parent)))<obj.limitTime);
                at = at(good_idx,:);
            end
            
            obj.arrivalTable =at;
        end
        
        function drawAgents(obj)
            
            % Update the arrival array
            if isempty(obj.arrivalArray)
                UpdateArrArray(obj);
            end
            
            temp_array = obj.arrivalArray;
            temp_array(:,end) =obj.Cluster_id;
            % Get rid of rows without a location
            good_idx = logical(~isnan(temp_array(:,end-1)));
            temp_array = temp_array(good_idx,:);
            
            
            duration = temp_array(end,1)-temp_array(1,1);
            newTimescale = round(temp_array(:,1)-temp_array(1,1))+1;
            
            
            
            % Hydrophone locations
            hyd_table = struct2table(obj.hydrophone_struct);
            TimeColorVals = parula(duration+5);
            ColorVals = lines(max(temp_array(:,end)));
            
            child_hyds = obj.array_struct.slave(obj.child_idx)
            
            
            figure;
            subplot(2,1,1)
            hold on
            scatter(temp_array(:,end-1),...
                temp_array(:,end-2),[],...
                TimeColorVals(newTimescale,:), 'f')
            scatter(hyd_table.location(:,2),...
                hyd_table.location(:,1), 80,...
                'k', 'filled', 'd')
            scatter(...
                hyd_table.location([obj.array_struct.master,...
                obj.array_struct.slave(obj.child_idx)],2),...
                hyd_table.location([obj.array_struct.master,...
                obj.array_struct.slave(obj.child_idx)],1),...
                'r', 'filled', 'd')
            
            title('TOA on Parent')
            
            subplot(2,1,2)
            hold on
            scatter(temp_array(:,end-1),...
                temp_array(:,end-2),[],...
                ColorVals(temp_array(:,end),:), 'f')
            scatter(hyd_table.location(:,2),...
                hyd_table.location(:,1), 80,...
                'k', 'filled', 'd')
            scatter(...
                hyd_table.location([obj.array_struct.master,...
                obj.array_struct.slave(obj.child_idx)],2),...
                hyd_table.location([obj.array_struct.master,...
                obj.array_struct.slave(obj.child_idx)],1),...
                'r', 'filled', 'd')
            
            
            
            
            
            
            
        end
    end
end