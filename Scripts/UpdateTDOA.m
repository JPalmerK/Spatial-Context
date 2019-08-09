        function TDOA_vals = UpdateTDOA(simStruct)
            % Function for extracting/calculating TDOA values from arrival
            % times matrix. In progress- allow different array configurations
            
            % Check if the arrival table is present if not update it
            if isempty(simStruct.arrivalArray)
                disp('Updating Arrival Array')
                UpdateArrArray(simStruct)
                
            end
            TDOA_vals =[];
            
            % Time difference of arrivals (can only handle two atm)
            
            for jj =1:length(simStruct.child_idx)
                TDOA_vals = [TDOA_vals,...
                    simStruct.arrivalArray(:, jj+1)-simStruct.arrivalArray(:, 1)];
            end
            
        end
        