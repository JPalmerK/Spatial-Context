        function array = UpdateArrArray(simStruct)
            % This funcion creates the array of arrivals, if the
            % arrivalsMiss parameter is set to 1, then where multiple calls
            % fall within the expected TDOA range, the call with the
            % arrival TDOA pair nearest to TDOA of 0 is chosen.
            % This selection is somewhat arbritrary as per design.
            
            % Check if the arrival table is present if not update it
            if isempty(simStruct.arrivalTable)
                disp('Updating Arrival Table')
                UpdateArrTable(simStruct)
                
            end
            
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
            
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Create the Associations and TDOA values %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % If clock drift is present, add/subtract it from the channels
            if isfield(simStruct, 'drift')
                if simStruct.drift ~=0
                    disp(['Add/Sub ', num2str(simStruct.drift), ' s of clock drift'])
                    for channum =1:length(simStruct.child_idx)
                        
                        % Shift the arrivals to simulate clock drift by either
                        % positive or negative values
                        clock_shift = (mod(channum, 2)*2-1) * simStruct.drift;
                        
                        array(:,channum+1) = array(:,channum+1) + clock_shift;
                    end
                end
            end
        
        
        
            
            % If we are adding in random missassociation then modify the
            % arrival array
            if isfield(simStruct, 'assSec') 
                
                disp(['Allowing for ', num2str(simStruct.assSec),...
                    ' s of misassociation' ])
                
                % Create a matrix that incorporates miss association
                ArrivalsMiss = zeros(size(array));
                ArrivalsMiss(:,1)  = array(:,1);
                ArrivalsMiss(:,end)  = array(:,end);
                
                % step through each child hydrophone and caluculate the
                % TDOA and maximum expected TDOA based on the hydrophone
                % array
                for channum =1:length(simStruct.child_idx)
                    
                    
                    % Allow for Mis Association %%%%%%%%%%%%%%%%%
                    
                    % Distance between two sensors - maximum exptected TDOA
                    % depth between calling whale and array
                    depth_range = simStruct.hydrophone_struct(simStruct.array_struct.master).depth...
                        - simStruct.hydrophone_struct(simStruct.array_struct.slave(channum)).depth ;
                    
                    % Calculate the horizontal distance between the calling
                    % locations and the hydrophone
                    horizontal_distance = arrayfun(@(lats, lons)...
                        vdist(lats, lons,...
                        simStruct.hydrophone_struct(simStruct.array_struct.master).location(1),...
                        simStruct.hydrophone_struct(simStruct.array_struct.master).location(2)),...
                        simStruct.hydrophone_struct(simStruct.array_struct.slave(channum)).location(1),...
                        simStruct.hydrophone_struct(simStruct.array_struct.slave(channum)).location(2));
                    
                    
                    % Maximum expected tdoa between parent and child
                    % (channum)
                    MaxTOA_1 = sqrt(depth_range^2 + ...
                        horizontal_distance.^2)/simStruct.c;
                    
                    
                    % Step through the arrivals and select a call at random
                    % from all calls within the MaxTOA_1
                    for ii = 1:length(ArrivalsMiss)
                        
                        % Get index of all calls that fall within the
                        % association zone plus the wiggle room
                        corr_idxs = find(...
                            abs(array(:,(channum+1))-array(ii,1))<...
                            MaxTOA_1+simStruct.assSec);
                        
                        % If more than one index falls in the correlation
                        % zone pick one
                        if length(corr_idxs)>1
                            
                            
                            % Random correlation
                            corr_idxs = datasample(corr_idxs,1);
                            
                            
                        elseif isempty(corr_idxs)
                            continue;
                        end
                        ArrivalsMiss(ii, channum+1) = array(corr_idxs,...
                            channum+1);
                    end
                    
                end
                
                % Change the 0's to nans to make future life easier
                ArrivalsMiss(ArrivalsMiss == 0) = nan;
                array = ArrivalsMiss;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
        end
        