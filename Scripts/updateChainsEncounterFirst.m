   function cluster = updateChainsEncounterFirst(obj)
            
            % Original design via Glen Lerley
            if isempty(obj.Sim_mat)
                disp('Updating Simulation Matrix')
                UpdateSimMat(obj)
            end
            
            

            
            % Ladder Linkages rebuild
            simmat= obj.Sim_mat;
            arraivalArray = obj.arrivalArray;
            max_gap = obj.maxEltTime;
            simthresh = obj.cutoff
            
            clusterN=1;
            
            cluster = struct();
            
            % Index of where we are in the array
            idx = 1:size(arraivalArray,1);
            
            % For each row, find the first gap
            while size(simmat,1)>1
                
                col_idx = idx(1);
                
                simscores = simmat(1,col_idx:end);
                times = arraivalArray(col_idx: col_idx+ sum(~isnan(simscores)));
                
                
                call_spacing = diff(times);
                
                % Clusters end at the first gap larger than maximum elapsed time
                can_clusterend = min([length(call_spacing), find(call_spacing>max_gap,1)]);
                
                % If the the first gap larger than the maximum elapsed time it gets
                % it's own cluster
                if can_clusterend<1
                    cluster(clusterN).index = col_idx;
                    cluster(clusterN).n =1;
                    
                    % Remove the rows from thematrix 
                    simmat(1,:)=[];
                    idx(1)=[];
                    
                    
                    
                else
                    % otherwise get the scores and remove the ones that are above the threshold
                    
                    matching_idx = unique([1, find(simscores(1:can_clusterend)>=simthresh)]);
                    
                    % If it only matches with it self, make a new cluster
                    if isempty(matching_idx)
                        
                        cluster(clusterN).index = col_idx;
                        
                        cluster(clusterN).n =1;
                        
                        % Remove the rows from thematrix
                        simmat(1,:)=[];
                        idx(1)=[];
                        
                    else
                        cluster(clusterN).index = col_idx+matching_idx-1;
                        cluster(clusterN).n = length(matching_idx);
                        disp(col_idx)
                        

                        % Remove the rows from thematrix
                        simmat(matching_idx,:)=[];
                        simmat(:,matching_idx)=[];
                        idx(matching_idx)=[];
                    end
                    
                    
                end
                
                clusterN=clusterN+1;
                
            end
            
            
 
            % aggregate clusters
            %chain = cluster;
            cluster
        end
        