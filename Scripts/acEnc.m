        %% Cluster based only on time of arrivals Baseline (step 4) equivallent to acoustic encounters
        function cluster_vals = acEnc(obj)
            
            % Check if the arrival table is present if not update it
            if isempty(obj.arrivalArray)
                disp(['Updating arrival Array'])
                UpdateArrArray(obj);
                
            end
            % Look for gaps bigger than the maximum overlap time
            diff_vals = diff(obj.arrivalArray(:,1));
            
            %             % Find gaps bigger than the allowable time
            %             time_idx = quantile(diff_vals, .95);
            %
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            % For sensitivity analysis
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            time_idx = obj.maxEltTime;
            
            idx = [1; find(diff_vals>time_idx)];
            
            
            cluster_vals= zeros(size(obj.arrivalArray(:,1)))+length(idx);
            cluster_id =1;
            
            
            for ii = 2:length(idx)-1
                
                cluster_vals(idx(ii-1):idx(ii)) = cluster_id;
                cluster_id = cluster_id +1;
                
            end
            
            
            cluster_vals;
            obj.titleStr = 'Baseline - toa Only';
            
        end