 function simMat= simMatIdealXcorrDist(simStruct)
            % This function creates the simulation matrix using the low
            % memory approach. This should be used in most cases where not
            % exploring the algorithims in depth.
            
            % Check if the arrival table is present if not update it
            if isempty(simStruct.TDOA_vals)
                disp(['Updating TDOA values'])
                UpdateTDOA(simStruct);
                
            end
            
            simStruct.titleStr ='Call Space Similarity Ideal';
            
            % Grid X/Y space
            % Get distance in meters between the lower and upper right
            grid_v = vdist(min(simStruct.array_struct.latgrid),...
                min(simStruct.array_struct.longrid),...
                max(simStruct.array_struct.latgrid),...
                min(simStruct.array_struct.longrid));
            
            
            % Get distance in meters between the lower left and lower right
            grid_h = vdist(min(simStruct.array_struct.latgrid),...
                min(simStruct.array_struct.longrid),...
                min(simStruct.array_struct.latgrid),...
                max(simStruct.array_struct.longrid));
            
            
            % Create the empty similarity matrix
            Sim_mat = zeros(length(simStruct.TDOA_vals))/0;
            
            
            % Grid X/Y space
            deltalat_space = (grid_v/ (length(simStruct.array_struct.latgrid)-1));
            deltalon_space = (grid_h/ (length(simStruct.array_struct.longrid)-1));
            
            % How many grid squares per second can the whale move
            lat_persec = simStruct.s / deltalat_space;
            lon_persec = simStruct.s / deltalon_space;
            
            
            sig_tot = sqrt(simStruct.PosUncertsigma + simStruct.drift^2);
            %profile on
            % Step through each arrival and get it's grid probability as
            % well as the projected grid probabilities for times at all
            % subsiquent calls but within the maximum time cuttoff
            simStruct.arrivalArray= (simStruct.arrivalArray);
            
            
                
           
            for ii =1:size(simStruct.arrivalArray,1)
                
                % Get the average prob loc space of the i-th call with
                % delta sigma t
                
                averageLklhd_space = (getTruHdSpace(simStruct, ii, sig_tot));
               
                % Figure out the number of time gaps within the maximum
                % allowed correlation time (time_cut)
                time_gaps = simStruct.arrivalArray(ii:end, 1)-...
                    simStruct.arrivalArray(ii, 1);
                
                diff_times = diff(simStruct.arrivalArray(ii:end, 1));
                % Find first big gap
                idx_end = find(diff_times>= simStruct.maxEltTime,1);
                
                time_gaps = time_gaps(1:idx_end);
                % If there are more than one time gap over which we need to
                % look then do the projections

                if length(time_gaps)>1
                        
                    
                    % Step through the time gaps/sigma values getting each
                    % probability loc space and projection
                    simValue = [];
                    
                    
                        % Grow the likelihood space based using image
                        % processing max filter. Set the filter size based
                        % on the maximum swim speed
                        filt_size_lat = round(lat_persec * time_gaps);
                        filt_size_lon = round(lon_persec * time_gaps);
                        filt_size = [filt_size_lat'; filt_size_lon'];
                        
                        Lklhd_space_proj_out = ones([...
                            length(simStruct.array_struct.latgrid),...
                            length(simStruct.array_struct.longrid),...
                            length(filt_size)]);
                        Lklhd_space_proj_out(:,:,1) = averageLklhd_space;
                        
   
                        for kk=2:length(filt_size)     
                            Lklhd_space_proj_out(:,:,kk) = imdilate(averageLklhd_space,...
                            true(filt_size(:,kk)'));
                            
                        end
                        

                     
                    for jj= 1:length(time_gaps)
                        
                        % If there a filter then project the space,
                        % otherwise don't. use 3d max filter based on the
                        % time gaps

                        Lklhd_space_proj=(...
                            squeeze(Lklhd_space_proj_out(:,:,jj)));
                        
                        % Get the prob. loc space space of the next call in
                        % the series
                        nextLklhdSpace = ...
                            (getTruHdSpace(simStruct,(ii+jj-1), sig_tot));
                        
                        aa = ((nextLklhdSpace));
                        bb = ((squeeze(Lklhd_space_proj_out(:,:,jj))));
                        latp = (gather(deltalat_space));
                        lonp = (gather(deltalon_space));
                        
                        %Stic
                        [dist, sim] = crossCorrSimScores(simStruct,aa,bb,...
                            latp, lonp);
                        speed =dist/time_gaps(jj);
                        %toc
                        
                        if speed>simStruct.s
                            sim=0;
                        end
                        
                        
                        simValue= [simValue sim];
                       
                        
                    end
                    
                    
                    
                    Sim_mat(ii, ii:ii+length(simValue)-1)=simValue;
                    Sim_mat(ii:ii+length(simValue)-1,ii)=simValue;
                    
                disp([num2str(ii), ' of ', num2str(length(simStruct.arrivalArray))])
                end
                %
%                 % Populate the simulation matrix
%                 Sim_mat(ii, ii:ii+length(simValue)-1) = simValue;
%                 Sim_mat(ii:ii+length(simValue)-1,ii) = simValue;
                %profile report
                %profile off
                %disp([num2str(ii), ' of ',...
                %    num2str(length(obj.arrivalArray))])
            end
            
          simMat  =Sim_mat;
            
        end