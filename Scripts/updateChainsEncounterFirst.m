   function cluster = updateChainsEncounterFirst(obj)
            
            % Original design via Glen Lerley
            if isempty(obj.Sim_mat)
                disp('Updating Simulation Matrix')
                UpdateSimMat(obj);
            end
            
            

            
            % Ladder Linkages rebuild
            simmat= obj.Sim_mat;
            arraivalArray = gather(obj.arrivalArray);
            max_gap = obj.maxEltTime;
            simthresh = obj.cutoff;
            
            clusterN=1;
            
            cluster = struct();
            
            % Index of where we are in the array
            idx = 1:size(arraivalArray,1);
            
            % For each row, find the first gap
            while size(simmat,1)>1
                
                col_idx = idx(1);
                
                simscores = simmat(1,:);
                % If the the first gap larger than the maximum elapsed time it gets
                % it's own cluster
                
                % otherwise get the scores and remove the ones that are above the threshold
                
                matching_idx = unique([1, find(simscores>=simthresh)]);
                
                
                cluster(clusterN).index = col_idx+matching_idx-1;
                cluster(clusterN).n = length(matching_idx);
                
                
                
                % Remove the rows from thematrix
                simmat(matching_idx,:)=[];
                simmat(:,matching_idx)=[];
                idx(matching_idx)=[];

                
                clusterN=clusterN+1;
            end
            
                
                
            end
            
            
 

       
        