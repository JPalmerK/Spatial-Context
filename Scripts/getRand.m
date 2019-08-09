function AdjRand= getRand(obj)
            
            % Check if cluster id's are present if not update it
            if isempty(obj.Cluster_id)
                updateClusterID(obj)
            end
            
            
            % Pull out the arrival array for easy reading
            arrival_array = obj.arrivalArray;
            trueClusters =unique(arrival_array(:,end));
            
            % Create smaller clusters based on the cut time
            timeClusters = acEnc(obj);
            
            
            % Create a dummy varaible for the clusters
            NewTruth = (arrival_array(:,end));
            
            
            % step through each true cluster and determine if it crosses a
            % break. If so it gets a new name;
                nextIdx = max(trueClusters)+1;
                for trueClust = 1:length(trueClusters)
                    
                    % Indicies of the true cluster
                    idx = find(arrival_array(:,end)== trueClust);
                    
                    % Predicted values based on time only
                    timeClustsub = timeClusters(idx);
                    
                    % List of time clusters
                    timClustList = unique(timeClustsub);
                    
                    % If the true cluster is in two temporal clusters, then
                    % adjust the cluster indexes
                    if length(timClustList)>1
                        
                        for jj=2: length(timClustList)
                            
                            % Index of the origional calls
                            new_idx = idx(find(timeClustsub== timClustList(jj)));
                            NewTruth(new_idx) = nextIdx;
                            nextIdx =nextIdx+1;
                            
                        end
                        
                    end
                    
                end

            
            
            
            
            
            % If the arrival array isn't empty, grab the adjusted rand idx
            if ~isempty(arrival_array)
                
                % Allign the true clusters and the predicted clusters
                newClusterId = allignclusters(NewTruth,...
                    obj.Cluster_id)+1;
                
                % Get the adjusted rand (third party)
                [f,~,~,~] = RandIndex(newClusterId, NewTruth);
                AdjRand = f;
            else
                % Otherwise, there were too few data to cluster
                AdjRand = nan;
                disp('No rand values available')
            end
            
            
        end