%% Spatial Context run system class %%
% 
%%

classdef clusterDataGPLclass <handle
    
    properties
        time_cut = 15*60 % seconds
        truncateKm = 15 % Distance beyond which calls are not detected
        hydrophone_struct % Hydrophone structure
        Parent_hyd = 19;
        CrossScore = .8
        array_struct % Array structure
        localize_struct % Localized calls from GPL system
        NChidlHyd = 2 % Minimum number of hydrophones over which to look at range
        s = 12 % Maximum swim speed (m/s)
        c = 1500 % speed of sound (m/s)
        fs =2000;
        cutoff = 0.6; % Correlation threshold cut off
        child_idx =[2,3]; % Column index for the child arrays
        high_prob_threshold =0.8; % High probability threshold for comparing two LSQ spaces
        arrivalTable
        arrivalArray
        TDOA_vals
        projSpace
        Sim_mat
        chains
        Cluster_id
        AdjRand
        truthAndpreds
        
    end
    
    
    methods
                %% Function for calculating the averageed (across hydrophones) 
        % probloc for a given call
        
        function averageLklhd_space = getAvLkHdSpace(obj, callIdx)
            % Returns an averaged likelihood space for a given TDOA
            
            % Pick the set of call delays (number of calls)
            delays = (obj.TDOA_vals(callIdx,:));
            
            % Remove nans, but keep the indexes for hydrophone pairs
            delays = delays(~isnan(delays));
            
            averageLklhd_space = [];
            
            % Iterate through the hydrophone pairs and get the combined TDOA
            for jj = 1:length(delays)
                
                % Get the tdoa space
                toa_space = (cell2mat(obj.array_struct(1).toa_diff(jj+1)));
                
                % potential toa space
                posLocSpace = (toa_space - delays(jj));
                
                % Convert tdoa space to probability space
                averageLklhd_space(:,:,jj) = posLocSpace;
                
            end
            
            % If there were multiple arrivals of the same call on each hydrophone,
            % take the average of the normalized PDF space
            if length(delays)>1
                
                % sum along third axis, will be normalized later
                averageLklhd_space = mean(averageLklhd_space,3);
                
            end
            
        end
        %% Create all habitat/area projections within the time cut (step 4)
        function simMatLowMemory(obj)
                      % This function creates the simulation matrix using the low
            % memory approach. This should be used in most cases where not
            % exploring the algorithims in depth. 
            
            % Check if the arrival table is present if not update it
            if isempty(obj.TDOA_vals)
                disp(['Updating TDOA values'])
                UpdateTDOA(obj);
                
            end
            
            % Create the empty similarity matrix
            Sim_mat = zeros(length(obj.TDOA_vals))/0;
            
            % Step through each arrival and get it's grid probability as
            % well as the projected grid probabilities for times at all
            % subsiquent calls but within the maximum time cuttoff
            for ii =1:length(obj.arrivalArray)
                Sim_mat(ii ,ii ) =1;
                
                % Get the average prob loc space of the i-th call 
                averageLklhd_space = getAvLkHdSpace(obj, ii); % need to rename this fx...
                
                % Figure out the number of time gaps within the maximum
                % allowed correlation time (time_cut)
                time_gaps = obj.arrivalArray(ii:end, 1)-...
                    obj.arrivalArray(ii, 1);
                time_gaps = time_gaps(time_gaps<obj.time_cut);
                
                % If there are more than one time gap over which we need to
                % look then do the projections
                
                if length(time_gaps)>1
                    
                    % Calculate the sigma values based on EM Nosal
                    % suggestion
                    % sigma (error) valuse from the normal distribution
                    sigma = 0.5+(time_gaps(2:end) * (2*obj.s)/obj.c);
                    
                    % Step through the time gaps/sigma values getting each
                    % probability loc space and projection
                    for jj= 1:(length(sigma)-1)
                        
                        % Grow the prob loc space of the current call
                        Lklhd_space_proj =  normpdf((averageLklhd_space),...
                            0, sigma(jj));
                        
                        % Get the prob. loc space space of the next call in
                        % the series 
                        nextLklhdSpace = getAvLkHdSpace(obj, (ii+jj));
                        nextLklhdSpace =  normpdf(nextLklhdSpace, 0, 0.5+(2*obj.s)/obj.c);
                        
                        % Get the comparison value of the projected prob
                        % loc space and the next call in the squence
                        simValue = compareLklhdSpace(obj, Lklhd_space_proj,...
                            nextLklhdSpace);
                        
                        % Populate the simulation matrix 
                        Sim_mat(ii, ii+jj ) = simValue;
                        Sim_mat(ii+jj ,ii) = simValue;
                        
                    end
                    
                end
                %disp([num2str(ii), ' of ',...
                %    num2str(length(obj.arrivalArray))])
            end
            obj.Sim_mat= Sim_mat;
            
            % Plot of the similarity matrix for good meausre
            figure;
            mycolormap = [ ones(1,3); (jet(30))];
            imagesc(Sim_mat)
            colormap(mycolormap);
            colorbar
            title('Call Space Similarity')
            xlabel('Call ID')
            ylabel('Call ID')
        end
        %% Create function to get TDOA values (step 2)
        function UpdateTDOA(obj)
            % Function for extracting/calculating TDOA values from arrival
            % times matrix. TO DO- allow different array configurations
            
            % Check if the arrival table is present if not update it
            if isempty(obj.arrivalArray)
                disp('Updating Arrival Array')
                UpdateArrArray(obj)
                
            end
            
            % Time difference of arrivals (can only handle two atm)
            TDOA_vals = [obj.arrivalArray(:, 1)- ...
                obj.arrivalArray(:, obj.child_idx(1)),...
                obj.arrivalArray(:, 1)-...
                obj.arrivalArray(:, obj.child_idx(2))];
            
            obj.TDOA_vals = TDOA_vals;
        end
        
        %% Create array of arrivals with or w/o random association (step 2) 
        function UpdateArrArray(obj)
            % This funcion creats the array of arrivals, if the
            % arrivalsMiss parameter is set to 1, then where multiple calls
            % fall within the expected TDOA range-, the call with the 
            % arrival TDOA pair nearest to TDOA of 0 is chosent.
            % This selection is somewhat arbritrary as per design.
            
            % Check if the arrival table is present if not update it
            if isempty(obj.arrivalTable)
                disp('Updating Arrival Table')
                UpdateArrTable(obj)
                
            end
            
            
            % Array containing the parent hyd (first column)and two selected child hyd
            array = [obj.arrivalTable.ArrivalSec(:,[1 obj.child_idx])];
            
            % Remove calls where there there are insufficient number of arrivals
            array = array(sum(~isnan(array(:,1:3)),2)>obj.NChidlHyd,:);
            
            % Sort by start time
            [~,sortidx] = sort(array(:,1));
            array = array(sortidx,:);
            obj.arrivalArray =array;
        end                
        
        %% Create table of arrival times (step 1)
        function UpdateArrTable(obj)
            % Use the TDOA table to create a matrix with TDOAs that are
            % clipped based on cross correlation score
 
            
            ParentArrival = obj.localize_struct.hyd(obj.Parent_hyd).rtimes./obj.fs;
            ParentArrival = [ParentArrival-min(ParentArrival)]';
            
            delays =obj.localize_struct.hyd(obj.Parent_hyd).delays;
            
            ChildArrival  = repmat(ParentArrival, 1,size(obj.localize_struct.hyd(obj.Parent_hyd).delays, 2));
            ArrivalSec = [ParentArrival ((delays+ChildArrival))];
            
                x=squeeze(obj.localize_struct.hyd(obj.Parent_hyd).coordinates(2,2,:));
                x(:,2)=squeeze(obj.localize_struct.hyd(obj.Parent_hyd).coordinates(2,1,:));
            
            % Create the arrivals table (at) from the localize structure
            at = table(...
                ArrivalSec,...
                x,...
                obj.localize_struct.hyd(obj.Parent_hyd).cross_score,...
                delays,...
                'VariableNames',{'ArrivalSec', 'Loc ','CrossScore', 'delays'});
            

            
            % Clear out any detections greater than Xkm from the receiver
            crossScore_bool = ones(size(at.CrossScore));
            crossScore_bool(find(at.CrossScore>obj.CrossScore)) = nan;
            
            % Add an extra column for the parent hydropone
    
            
            at.CrossScore = at.CrossScore.*crossScore_bool;
            at.ArrivalSec(:,2:end) = at.ArrivalSec(:,2:end).*crossScore_bool;
            obj.arrivalTable =at;
        end
        
        %% Clear calculated values (simmat, projections, TODA values etc)
        function clearCalcValues(obj)
            
            % Clear out all values that need to be calculated leaving the
            % parameters intact. Useful for running multiple experiments
            % using the same agent-setup but different correlation or gap
            % parameters
                obj.arrivalTable=[];
                obj.TDOA_vals=[];
                obj.projSpace=[];
                obj.Cluster_id=[];
                obj.AdjRand=[];
                obj.Sim_mat =[];
                obj.arrivalArray =[];
                obj.chains =[];
                obj.AdjRand=[];
                obj.truthAndpreds =[];
             

        end
        
        %% Compare two probability grid spaces (Help me!)
        function simValue = compareLklhdSpace(obj,...
                ProjectedLklhdSpace, nextLklhdSpace)
            
            % Compares likelihood spaces for two calls, typically call a that
            % has been projected and call b from later in the sequence that
            % has not
           
            % Inputs - 
            % ProjectedLklhdSpace - likelihood map for the call projected
            % across time
            % nextLklhdSpace - likelihood map for the next call in the
            % sequence (typically)
            % Returns-
            % 
            
            % Normalize the space
            nextLklhdSpacenorm = (nextLklhdSpace - min(min(nextLklhdSpace)))./...
                ( max(max(nextLklhdSpace)) - min(min(nextLklhdSpace)));
            
            % Normalize the space
            ProjectedLklhdSpace = (ProjectedLklhdSpace - ...
                min(min(ProjectedLklhdSpace)))./...
                ( max(max(ProjectedLklhdSpace)) -...
                min(min(ProjectedLklhdSpace)));
            
%             
%             % Size of area represented by the next call
%             b_hiprob = find(nextLklhdSpace> obj.high_prob_threshold);
%             size_b_highprob = length(b_hiprob);
%             
%             a_hiprob = find(ProjectedLklhdSpace> obj.high_prob_threshold);
%             
%             
%             % Determine the percentage of the next call that is covered by
%             % the time expanded current call 
%             simValue = length(intersect(a_hiprob, b_hiprob))/size_b_highprob;
%             
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            simValue = sum(sum(nextLklhdSpacenorm.*ProjectedLklhdSpace))/...
                sum(sum(nextLklhdSpacenorm));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
        end

        %% Function for creating predicted cluster ID's for the agents
        function updateClusterID(obj)
            
            % If chains have not been updated, do so
            if isempty(obj.chains)
                disp('Updating Simulation Matrix')
                updateChains(obj)
            end
            
            % Simply create an array with the predicted clusterid
            Cluster_id = zeros(length(obj.TDOA_vals),1);
            
            for ii=1:length(obj.chains);
                Cluster_id(obj.chains(ii).index) = ii;
                
            end
            
            obj.Cluster_id =Cluster_id;
        end
        
        %% Function for creating the cluster chains (aka ladder linkages)
        function updateChains(obj)
            
            if isempty(obj.Sim_mat)
                disp('Updating Simulation Matrix')
                UpdateSimMat(obj)
            end
            
            
            % Ladder linkages for new reduced plot
            cutoff =obj.cutoff; % Correlation threshold cut off
            time_cut = obj.time_cut;
            fs = obj.fs
            

            time = obj.localize_struct.hyd(obj.Parent_hyd).rtimes; 


            % Elapsed time since first call
            time=(time-time(1))/fs;  % seconds now yes?
            time0=time; % keep original copy for later indexing


            time = obj.arrivalArray;
            
            % Elapsed time since first call
            time=(time-time(1));  % seconds
            time0=time; % keep original copy for later indexing
            
            % Replicate Corr_coef_map for shrinking matrix
            Corr_coef_map = obj.Sim_mat;
            [s1,s2]=size(Corr_coef_map);
            nc=1; % cluster counter
            clear chain dex ss
            
            while s1>1
                flag=1;
                
                % call element index, start at 1 and examine columns
                elt=1; % first element of cluster is always upper left hand element
                index=[];
                cluster=[];
                
                
                % Find all calls that are sufficiently correlated with this one.
                while flag == 1
                    % Add call to current cluster, retaining the time and index
                    % of the call element
                    cluster = [cluster, time(elt)];
                    index   = [index,elt];
                    
                    % entire column for each element
                    vec=Corr_coef_map(:,elt);
                    
                    % Only look at the elements below
                    vec=vec(elt+1:end);
                    
                    % Look for correlations that exceed the threshold
                    k=find(vec>cutoff);
                    
                    if isempty(k)
                        flag=0;  % Nobody met the criterion
                    else
                        elt=k(1)+elt;  % Set next call to the first that meets criterion
                    end
                end % flag test
                
                % now cut cluster at first large jump:
                k=find(diff(cluster)>time_cut);
                if length(k)>0
                    k=k(1);
                    % Remove everything after the first large gap
                    cluster=cluster(1:k);
                    index=index(1:k);
                end
                
                % aggregate clusters
                chain(nc).links=cluster;
                chain(nc).n = length(cluster);
                for j=1:length(cluster)
                    
                    cluster_j = min(find(time0==cluster(j)))
                    chain(nc).index(j)=cluster_j;
                
                end
                nc=nc+1;
                
                % filter out this cluster to reduce matrix
                [s1,s2]=size(Corr_coef_map);
                rs1=setdiff([1:s1],index);
                Corr_coef_map=Corr_coef_map(rs1,:);
                Corr_coef_map=Corr_coef_map(:,rs1);
                time=time(rs1);
                [s1,s2]=size(Corr_coef_map);
                
            end % while s1 loop
            
            chain = struct(chain);
            
            obj.chains =chain;
        end
        
        %% Create all habitat/area projections within time cut (RAM heavy!)
        function UpdateprojSpace(obj)
            % This is a previous version of the simMatLowMemory function. 
            % It saves all projected spaces in the array and as such often 
            % causes the program to crash due to memory requirements. 
            % However, it is useful for saving and investigating differen
            % projected and non/habitat spaces for making figures and 
            % Developing comparison algorithithms. 
            
            % Check if the arrival table is present if not update it
            if isempty(obj.TDOA_vals)
                disp(['Updating TDOA values'])
                UpdateTDOA(obj);
                
            end
            
            % Create structure for each call and it's projected space
            Proj_space = struct();
            
            % Step through each arrival and get it's grid probability as
            % well as the projected grid probabilities for times at all
            % subsiquent calls but within the maximum time cuttoff
            for ii =1:length(obj.arrivalArray)
                
                % Pick a set of call delays (number of calls)
                delays = (obj.TDOA_vals(ii,:));
                
                % Remove nans
                delays = delays(~isnan(delays));
                
                averageLSQ_space = [];
                
                % Iterate through the hydrophone pairs and get the combined TDOA
                for jj = 1:length(delays)
                    
                    % Get the tdoa space
                    toa_space = (cell2mat(obj.array_struct(1).toa_diff(jj+1)));
                    
                    % potential toa space
                    LSQ_VAL = (toa_space - delays(jj)).^2;
                    
                    % Convert tdoa space to probability space
                    averageLSQ_space(:,:,jj) = LSQ_VAL;
                    
                end
                
                % If there were multiple arrivals of the same call on each hydrophone,
                % take the average of the normalized PDF space
                if length(delays)>1
                    
                    % average along the third axis
                    averageLSQ_space = mean(averageLSQ_space,3);
                    
                end
                
                % Do the time projection up to the maximum correlation time
                % figure out which calls are within the alotted fimeframe
                time_gaps = obj.arrivalArray(ii:end, 1)-...
                    obj.arrivalArray(ii, 1);
                
                % Only interested in time gaps less than the time cut off (5min)
                time_gaps = time_gaps(time_gaps<obj.time_cut);
                
                % sigma (error) valuse from the normal distribution
                sigma = (time_gaps * (2*obj.s)/obj.c);
                
                % Iterate through  the time gaps and create the projections
                for jj= 1:length(time_gaps)
                    
                    % Grow the LSQ space
                    LSQ_space =  normpdf(averageLSQ_space, 0, (1+sigma(jj)));
                    
                    Proj_space(ii).projection(:,:,jj) = LSQ_space;
                end
                if mod(ii,20)==0
                    %disp([num2str(ii), ' of ', num2str(length(obj.arrivalArray))])
                end
            end
            
            % Set the projection spaces property
            obj.projSpace = Proj_space;
        end
                        
    end
end

