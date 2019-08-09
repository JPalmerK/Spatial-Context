        %% Create table of arrival times (step 1)
        function at = UpdateArrTable(simStruct)
            % Use the TDOA table to create a matrix with TDOAs that are
            % clipped based on distance to receiver
            
            % Make life a little easier
            spaceWhaleloc = simStruct.spaceWhale;
            
            % Create the arrivals table (at) from the first agent
            at = table(...
                spaceWhaleloc.agent(1).call_times', ...
                spaceWhaleloc.agent(1).Arrival_times,...
                spaceWhaleloc.agent(1).RangeKm,...
                spaceWhaleloc.agent(1).location(spaceWhaleloc.agent(1).call_times,:),...
                ones(length(spaceWhaleloc.agent(1).call_times),1),...
                struct2array( spaceWhaleloc.agent(1).tdoa(simStruct.array_struct.master)),...
                'VariableNames',{'callTimes','ArrivalSec', 'RangeKm',...
                'Location','ID', 'TDOA'});
            
            % Fill in the arrivals table with data from the rest of the
            % agents
            for ii = 2:length(spaceWhaleloc.agent)
                Tnew =table(...
                    spaceWhaleloc.agent(ii).call_times', ...
                    spaceWhaleloc.agent(ii).Arrival_times,...
                    spaceWhaleloc.agent(ii).RangeKm,...
                    spaceWhaleloc.agent(ii).location(...
                    spaceWhaleloc.agent(ii).call_times,:),...
                    ones(length(spaceWhaleloc.agent(ii).call_times),1)*ii,...
                    struct2array( spaceWhaleloc.agent(ii).tdoa(simStruct.array_struct.master)),...
                    'VariableNames',{'callTimes','ArrivalSec', 'RangeKm',...
                    'Location','ID', 'TDOA'});
                at = [at ; Tnew];
            end
            
            % annotate out any detections greater than Xkm from the receiver
            distance_bool = ones(size(at.RangeKm));
            distance_bool(find(at.RangeKm>simStruct.truncateKm)) = nan;
            
            at.RangeKm = at.RangeKm.*distance_bool;
            at.ArrivalSec = at.ArrivalSec.*distance_bool;

        end
        