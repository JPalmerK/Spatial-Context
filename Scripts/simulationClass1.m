%% Spatial Context simulation class %%
% The following class allows users to simulate and measure the algorithism
% clustering ability based on a call detected on two or more sensors
% The code is desgned to work with the outputs of the GPL detector and the
% whales in space (space whales!) agent based movement models.
% The system is designed to ouput adjusted rand indices for the simmulated
% agents. Also, a simulated cluster decision classifier (Palmer et al 2017)
% is in cluded and the system will output the proportion of correct and
% misclassified calls.
% At present the system has a built in function to re-run the
%%

classdef simulationClass <handle
    
    % The data section
    properties
        time_cut = 10*60 % seconds - time beyond which not to corrlelations
        spaceWhale % Structure containing call times and locations of all agents
        truncateKm = 10 % Distance beyond which calls are excluded (not detected)
        hydrophone_struct % Hydrophone structure same format as GPL
        array_struct % Array structure - containing hydrophone pairs and
        % expected TDOA grids. GPL and localization format
        
        NChidlHyd = 2 % Minimum number of hydrophones over which to
        % look at. This is currently fixed at 2, but should be expanded
        % for one hydrophone and 2+ hydrophone pairs
        child_idx ; % Column index for the child arrays
        
        randomMiss =0 % Bool for random misassosiation (0) for ideal assoc.
        assSec = 4 % Number of seconds beyond the expected arraival time
        %over which to look for potentially missassociated calls.
        %Used only for creating arrival arrays with misassociation
        drift = 0 % Number of seconds to add/subtract to arrivals to
        % simulate clock drift
        
        s = 12 % Maximum swim speed of agent (m/s)
        c = 1500 % speed of sound (m/s)
        fs =2000; % sample rate
        
        %cutoff = 0.8; % Correlation threshold cutoff - how similar do call
        % spaces need to be in order to be clustered
        cutoff 
        maxEltTime  % maximum time between calls to start another cluster
        
        % The following variables are created as the system is run
        arrivalTable % Table containing arrival times, locations etc of
        % each call
        arrivalArray % Arrival times and ID's of each call
        TDOA_vals % TDOA values to feed to the clustering algorithim
        
        Sim_mat % Similarity matrix for the call probability spaces
        chains % Structure with all linked call indexes and number of calls
        % in each cluster
        Cluster_id % Array size of included agent calls indicating the
        
        % predicted cluster ID for each call
        tweakedRand =1% Bool for whether or not to adjust the truth values for the clusters.
        % Only reason you wouldn't do this is when looking for the true
        % rand.
        AdjRand % Adjusted rand index
        truthAndpreds % Matrix containing the call id's and the predicated
        % call ID's based on the algorithim
        
        wrongAssocIdx % Index of the calls with incorrect associations
        
        % Position/soundspeed uncertainty, see Eva's paper
        % receiver position = 5 m -> 5/1500 = 0.0004 sec
        % height LSQ peak (check with Tyler) = .1 sec
        % sound speed profile = 50 m/s -> 50/1500 =0.3 sec
        PosUncertsigma = 0.0004^2 +.1^2 + .3^2;
        
        % Logit of the correct classification rate for the
        betaParm1
        betaParm2
        
        
        % Title string for simulation plots
        titleStr
    end
    
    
    methods
        %% Assign score and species for each of the agents
        function createSpeciesPreds(obj)
            % This function is the classifier function that could use a bit
            % of love. It's a bit hacky. Looks at calls and assigns a
            % random probability score indicating that each call came from
            % species 1.
            
            % Check the arrival array, if it's empty populate it
            if isempty(obj.arrivalArray)
                UpdateArrArray(obj);
            end
            % Check to see if clusters id had been added
            if isempty(obj.Cluster_id)
                updateClusterID(obj);
            end
            
            truthAndpreds = obj.truthAndpreds;
            
            
            % Stick the real id's next to the predicted clusters
            truthAndpreds =[obj.arrivalArray(:,end) obj.Cluster_id];
            truthAndpreds(:,3)=0; % True class
            truthAndpreds(:,4)=0; % Score
            truthAndpreds(:,5)=0; % Prediction after applying likelihood ratio
            
            truthAndpreds= array2table(truthAndpreds, 'VariableNames',...
                {'TrueClust','PredClust', 'TrueSpp', 'Score',...
                'ClusterScore'});
            
            
            % Detector Parameters
            rw_mean = obj.betaParm1;
            rw_sd = obj.betaParm2;
            
            % For each detected agent, assign a species and classification
            % probability
            agent_idxs =unique(truthAndpreds.TrueClust);
            
            for ii =1:length(agent_idxs)
                
                % Agent id
                agent_id =agent_idxs(ii);
                
                % If it's odd it's definitely a humpback
                if mod(agent_id,2)
                    
                    % Get the indicies of the calls of that individual
                    rw_idx = find(truthAndpreds.TrueClust == agent_id);
                    
                    % Pull from distribution to estimate classificaiton probability
                    score = betarnd(obj.betaParm1, obj.betaParm2,[1 length(rw_idx)]);
                    %score = .6;
                    
                    
                    truthAndpreds.TrueSpp(rw_idx) =1; %yes target spp.!
                    truthAndpreds.Score(rw_idx) = score';
                    
                    
                else
                    
                    
                    % Same as above, only other species
                    
                    % Get the indicies of the calls of that individual
                    mn_idx = find(truthAndpreds.TrueClust == agent_id);
                    
                    % Transform so that it's between 0 and 1
                    score = 1- betarnd(obj.betaParm1, obj.betaParm2,[1 length(mn_idx)]); % logit (or inverse logit)
                    %score = .4;
                    
                    % Update the species (binary) and the classifier
                    % predication score
                    truthAndpreds.TrueSpp(mn_idx) = 0; % No! Not target spp.
                    truthAndpreds.Score(mn_idx) = score';
                    
                end
                
                
                
            end
            
            
            
            % Update the object with the new table
            obj.truthAndpreds =truthAndpreds;
            
        end
        
        %% Function for applying simulated detector/classifier
        function perfStructure = estClassifierPerf(obj)
            % Compare the classification performance of the system without
            % clustering to the system where all calls in a cluster are
            % classified together
            
            createSpeciesPreds(obj);
            
            
            % Pull out for easy reading
            truthAndpreds = obj.truthAndpreds;
            
            
            cluster_ids = unique(truthAndpreds.PredClust);
            % For each cluster create the classificcation balues
            for ii=1:length(cluster_ids)
                
                % get the ichain index value
                idx =find(truthAndpreds.PredClust==cluster_ids(ii));
                
                % Gt the classifier scores
                scores = truthAndpreds.Score(idx);
                
                
                % Likelihood ratio for the cluster
                LR = log(prod(scores./(1-scores)));
                
                % Use Geometric mean instead
                truthAndpreds.GeoMean(idx) = geomean(scores);
                
                % Fill in the likelihood ratio
                truthAndpreds.ClusterScore(idx) = LR;
                
                %Simple Voting
                pos_votes = sum(scores>0.5); neg_votes=sum(scores<=.5);
                
                if pos_votes>neg_votes
                    truthAndpreds.Voting(idx) = ones(size(scores));
                elseif pos_votes<neg_votes
                    truthAndpreds.Voting(idx) = zeros(size(scores));
                elseif pos_votes==neg_votes
                     truthAndpreds.Voting(idx) = zeros(size(scores))+binornd(1,.5);
                end
                
            end
            
            obj.truthAndpreds = truthAndpreds;
            perfStructure = table2struct(truthAndpreds);
            
            
        end
        
        %% Function for running the experiment with or/without random missclassification (steps 1-4)
        function runRandom(obj, val)
            % This function re-runs the simulation switching between ideal
            % and random association
            
            % If input value is the same as what it already is then check
            % whether it's been run. If not, run it.
            if val == obj.randomMiss
                
                % And chains is empty, run the experiment
                if isempty(obj.chains)
                    
                    % Clear the values that need to be recalculated
                    clearCalcValues(obj)
                    
                    % Otherwise no need to re-run
                else
                    disp(['Model has been run with random association set to ',...
                        num2str(obj.randomMiss)]);
                end
                
                % Otherwise, the model has been run and we need to re-run
                % based on the new input value
            else
                
                % Clear the values that need to be recalculated
                clearCalcValues(obj)
                
                % Update the object
                obj.randomMiss = val;
                
                
            end
            
            
        end
  %%      
        function simMatIdealNewSim(obj)
                      % This function creates the simulation matrix using the low
            % memory approach. This should be used in most cases where not
            % exploring the algorithims in depth.
            
            % Check if the arrival table is present if not update it
            if isempty(obj.TDOA_vals)
                disp(['Updating TDOA values'])
                UpdateTDOA(obj);
                
            end
            
            obj.titleStr ='Call Space Similarity Ideal';
            
            % Grid X/Y space
            % Get distance in meters between the lower and upper right
            grid_v = vdist(min(obj.array_struct.latgrid),...
                min(obj.array_struct.longrid),...
                max(obj.array_struct.latgrid),...
                min(obj.array_struct.longrid));
            
            
            % Get distance in meters between the lower left and lower right
            grid_h = vdist(min(obj.array_struct.latgrid),...
                min(obj.array_struct.longrid),...
                min(obj.array_struct.latgrid),...
                max(obj.array_struct.longrid));
            
            
            % Create the empty similarity matrix
            Sim_mat = zeros(length(obj.TDOA_vals))/0;
            
            
            % Grid X/Y space
            deltalat_space = grid_v/ (length(obj.array_struct.latgrid)-1);
            deltalon_space = grid_h/ (length(obj.array_struct.longrid)-1);
            
            % How many grid squares per second can the whale move
            lat_persec = obj.s / deltalat_space;
            lon_persec = obj.s / deltalon_space;
            
            sig_tot = 2*sqrt(obj.PosUncertsigma + obj.drift^2);
            
            % Step through each arrival and get it's grid probability as
            % well as the projected grid probabilities for times at all
            % subsiquent calls but within the maximum time cuttoff
            for ii =1:length(obj.arrivalArray)
                
                % Get the average prob loc space of the i-th call with
                % delta sigma t
                
                
                averageLklhd_space = getTruHdSpace(obj, ii, sig_tot);
                
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
                    
                    
                    % Step through the time gaps/sigma values getting each
                    % probability loc space and projection
                    for jj= 1:length(time_gaps)
                        
                       
                        nextLklhdSpace = getTruHdSpace(obj, (ii+jj-1), sig_tot);
                        
%                         % Get the comparison value of the projected prob
%                         % loc space and the next call in the squence
%                         aa = sum(averageLklhd_space,1);
%                         bb = sum(averageLklhd_space,2);
%                         
%                         
%                         aa2 = sum(nextLklhdSpace,1);
%                         bb2 = sum(nextLklhdSpace,2);
%                         
%                         cor1 = corrcoef(aa,aa2);
%                         cor2 = corrcoef(bb,bb2);
%                         simValue = (cor1(1,2) + cor2(1,2))/2;
                        [dist, corrScore]= crossCorrSimScores(...
                            averageLklhd_space,nextLklhdSpace,...
                            deltalat_space, deltalon_space);
                        
                        speed = dist/time_gaps(jj);
                        
                        % value is minimum of the normalized correlation
                        % threshold and maximum swim speed
                        
                        if speed<obj.s
                            simValue = corrScore;
                        else
                            simValue =0;
                        end
                        
                        
                        
%                         simValue = compareLklhdSpace(obj, Lklhd_space_proj,...
%                             nextLklhdSpace);
                        %                         figure
                        %                         subplot(2,1,1)
                        %                         imagesc(Lklhd_space_proj); colorbar
                        %
                        %                         subplot(2,1,2)
                        %                         imagesc(nextLklhdSpace); colorbar
                        %
                        
                        % Populate the simulation matrix
                        Sim_mat(ii, ii+jj-1) = simValue;
                        Sim_mat(ii+jj-1 ,ii) = simValue;
                        
                    end
                    
                end
                %disp([num2str(ii), ' of ',...
                %    num2str(length(obj.arrivalArray))])
            end
            obj.Sim_mat= Sim_mat;
            
            
            
        end
        %% Create all habitat/area projections within the time cut (step 4)
        function simMatIdeal(obj)
            % This function creates the simulation matrix using the low
            % memory approach. This should be used in most cases where not
            % exploring the algorithims in depth.
            
            % Check if the arrival table is present if not update it
            if isempty(obj.TDOA_vals)
                disp(['Updating TDOA values'])
                UpdateTDOA(obj);
                
            end
            
            obj.titleStr ='Call Space Similarity Ideal';
            
            % Grid X/Y space
            % Get distance in meters between the lower and upper right
            grid_v = vdist(min(obj.array_struct.latgrid),...
                min(obj.array_struct.longrid),...
                max(obj.array_struct.latgrid),...
                min(obj.array_struct.longrid));
            
            
            % Get distance in meters between the lower left and lower right
            grid_h = vdist(min(obj.array_struct.latgrid),...
                min(obj.array_struct.longrid),...
                min(obj.array_struct.latgrid),...
                max(obj.array_struct.longrid));
            
            
            % Create the empty similarity matrix
            Sim_mat = zeros(length(obj.TDOA_vals))/0;
            
            
            % Grid X/Y space
            deltalat_space = grid_v/ (length(obj.array_struct.latgrid)-1);
            deltalon_space = grid_h/ (length(obj.array_struct.longrid)-1);
            
            % How many grid squares per second can the whale move
            lat_persec = obj.s / deltalat_space;
            lon_persec = obj.s / deltalon_space;
            
            
            % Step through each arrival and get it's grid probability as
            % well as the projected grid probabilities for times at all
            % subsiquent calls but within the maximum time cuttoff
            for ii =1:length(obj.arrivalArray)
                
                % Get the average prob loc space of the i-th call with
                % delta sigma t
                
                sig_tot = 2*sqrt(obj.PosUncertsigma + obj.drift^2);
                averageLklhd_space = getTruHdSpace(obj, ii, sig_tot);
                
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
                    
                    
                    % Step through the time gaps/sigma values getting each
                    % probability loc space and projection
                    for jj= 1:length(time_gaps)
                        
                        
                        % Grow the likelihood space based using image
                        % processing max filter. Set the filter size based
                        % on the maximum swim speed
                        filt_size_lat = ceil(lat_persec * time_gaps(jj));
                        filt_size_lon = ceil(lon_persec * time_gaps(jj));
                        filt_size = max([filt_size_lat filt_size_lon]);
                        
                        % If there a filter then project the space,
                        % otherwise don't. use 3d max filter based on the
                        % time gaps
                        if filt_size>1
                            Lklhd_space_proj = imdilate(averageLklhd_space,...
                                true(filt_size));
                        else
                            Lklhd_space_proj =averageLklhd_space;
                        end
                        
                        % Get the prob. loc space space of the next call in
                        % the series
                        nextLklhdSpace = getTruHdSpace(obj, (ii+jj-1), sig_tot);
                        
                        % Get the comparison value of the projected prob
                        % loc space and the next call in the squence

                        
                        simValue = compareLklhdSpace(obj, Lklhd_space_proj,...
                            nextLklhdSpace);
                        %                         figure
                        %                         subplot(2,1,1)
                        %                         imagesc(Lklhd_space_proj); colorbar
                        %
                        %                         subplot(2,1,2)
                        %                         imagesc(nextLklhdSpace); colorbar
                        %
                        
                        % Populate the simulation matrix
                        Sim_mat(ii, ii+jj-1) = simValue;
                        Sim_mat(ii+jj-1 ,ii) = simValue;
                        
                    end
                    
                end
                %disp([num2str(ii), ' of ',...
                %    num2str(length(obj.arrivalArray))])
            end
            obj.Sim_mat= Sim_mat;
            
            
            
        end
        
        %% Create all habitat/area projections within the time cut (step 4)
        function simMatTDOAonly(obj)
            if isempty(obj.TDOA_vals)
                disp(['Updating TDOA values'])
                UpdateTDOA(obj);
            end
            obj.titleStr ='Call Space Similarity TDOA Only';
            % Create the empty similarity matrix
            Sim_mat = zeros(length(obj.arrivalArray(:,end)))/0;
            
            for ii =1:(length(obj.arrivalArray)-1)
                
                
                tdoa_orig = obj.TDOA_vals(ii,:);
                
                
                % index of good tdoa vals
                
                
                % index of all calls within the elapsed time
                nextTimes = obj.arrivalArray(ii:end,1);
                time_diffs =  nextTimes- obj.arrivalArray(ii,1);
                
                % index/trim diff times greater than the maximum elapsed time
                okDiffsIdx = find(time_diffs<= obj.time_cut);
                time_diffs = time_diffs(okDiffsIdx);
                
                % Get the TDOA values
                TDOA_next = obj.TDOA_vals((ii+okDiffsIdx-1),:);
                
                
                
                % For each hydrophone pair calculate likelihood values
                deltaTDOALklhd =[];
                
                
                for jj=1:size(TDOA_next,2)
                    
                    mu = zeros(length(time_diffs),1);
                    sigmaSwim = 2* sqrt((time_diffs * (obj.s)/obj.c).^2);
                    x =(tdoa_orig(jj) - TDOA_next(:,jj)); %values
                    
                    likelihood = normpdf(x,mu,sigmaSwim);
                    
                    % Normalizing factor
                    LikelihoodNormFac= normpdf(0,0,sigmaSwim);
                    NormLikelihood = likelihood./LikelihoodNormFac;
                    
                    
                    % Normalized likelihood
                    deltaTDOALklhd = [deltaTDOALklhd, NormLikelihood]; % sigma
                    % Create normalizing factor
                    
                    
                end
                
                simValues = nanmin(deltaTDOALklhd,[],2);
                
                % Take the minimum value and fill in the similarity matrix
                Sim_mat(ii,ii) = 1;
                Sim_mat(ii, ii:ii+length(simValues)-1) = simValues;
                Sim_mat(ii:ii+length(simValues)-1,ii) = simValues;
            end
            
            obj.Sim_mat= Sim_mat;
            
            
            
        end
        %% Create all habitat/area projections within the time cut (step 4)
        function simMatadHoc(obj)
            % This function creates the simulation matrix using the low
            % memory approach. This should be used in most cases where not
            % exploring the algorithims in depth.
            
            % Check if the arrival table is present if not update it
            if isempty(obj.TDOA_vals)
                disp(['Updating TDOA values'])
                UpdateTDOA(obj);
                
            end
            
            % Set title string for plots, also uesful in similarity matrix
            % comparion
            
            
            obj.titleStr = 'Call Space Similarity adHoc';
            
            % Create the empty similarity matrix
            Sim_mat = zeros(length(obj.arrivalArray))./0;
            
            
            % Step through each arrival and get it's grid probability as
            % well as the projected grid probabilities for times at all
            % subsiquent calls but within the maximum time cuttoff
            for ii =1:length(obj.arrivalArray)
                
                % Get the average prob loc space of the i-th call
                averageLklhd_space = getAvLkHdSpace(obj, ii); % need to rename this fx...
                
                %contourf(obj.array_struct.longrid, obj.array_struct.latgrid, averageLklhd_space1)
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
                    
                    swimSigma =  (time_gaps * ((obj.s)/obj.c));
                    
                    % Total sigma error (standard deviation) five input
                    % variables (four from position uncertainty plus sigma
                    % swim)
                    sigSD = 2*sqrt((obj.PosUncertsigma + swimSigma.^2+ obj.drift^2));
                    
                    % Step through the time gaps/sigma values getting each
                    % probability loc space and projection
                    for jj= 1:length(sigSD)
                        
                        % Grow the prob loc space of the current call and
                        % normalize
                        Lklhd_space_proj =  ...
                            normpdf(averageLklhd_space,0, sigSD(jj))./...
                            normpdf(0, 0, sigSD(jj));
                        
                        if ndims(Lklhd_space_proj)>2
                            Lklhd_space_proj= nanmean(Lklhd_space_proj,3);
                        end
                        
                        
                        % Get the prob. loc space space of the next call in
                        % the series and normalize
                        %
                        sigma = 2*sqrt(obj.PosUncertsigma+ obj.drift^2);
                        nextLklhdSpace = getTruHdSpace(obj, (ii+jj-1), sigma);
                        
                        
                        %
                        %                         close all; figure(1); subplot(2,1,1);
                        %                         imagesc(nextLklhdSpace); colorbar;
                        %                         subplot(2,1,2); imagesc(Lklhd_space_proj); colorbar
                        
                        % Get the comparison value of the projected prob
                        % loc space and the next call in the squence
                        simValue = compareLklhdSpace(obj, Lklhd_space_proj,...
                            nextLklhdSpace);
                        if isnan(simValue)
                            aa1
                        end
                        
                        % Populate the simulation matrix
                        Sim_mat(ii, ii+jj-1 ) = simValue;
                        Sim_mat(ii+jj-1 ,ii) = simValue;
                        
                    end
                    
                end
                %disp([num2str(ii), ' of ',...
                %    num2str(length(obj.arrivalArray))])
            end
            obj.Sim_mat= Sim_mat;
            
            
        end
        
        %% Cluster based only on time of arrivals Baseline (step 4)
        function cluster_vals = tempClust(obj)
            
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
        
        
        %% Baseline clustering against which all others are compared
        function toaOnlyCluster(obj)
            obj.Cluster_id = tempClust(obj);
        end
        
        %% Create function to get TDOA values (step 2)
        function UpdateTDOA(obj)
            % Function for extracting/calculating TDOA values from arrival
            % times matrix. In progress- allow different array configurations
            
            % Check if the arrival table is present if not update it
            if isempty(obj.arrivalArray)
                disp('Updating Arrival Array')
                UpdateArrArray(obj)
                
            end
            TDOA_vals =[];
            
            % Time difference of arrivals (can only handle two atm)
            for jj =1:length(obj.child_idx)
                TDOA_vals = [TDOA_vals,...
                    obj.arrivalArray(:, jj+1)-obj.arrivalArray(:, 1)];
            end
            obj.TDOA_vals = TDOA_vals;
        end
        
        %% Create array of arrivals with or w/o random association (step 2)
        function UpdateArrArray(obj)
            % This funcion creates the array of arrivals, if the
            % arrivalsMiss parameter is set to 1, then where multiple calls
            % fall within the expected TDOA range, the call with the
            % arrival TDOA pair nearest to TDOA of 0 is chosen.
            % This selection is somewhat arbritrary as per design.
            
            % Check if the arrival table is present if not update it
            if isempty(obj.arrivalTable)
                disp('Updating Arrival Table')
                UpdateArrTable(obj)
                
            end
            
            % Array containing the parent and children hydrophones
            
            % Get ID vlue, if simulated data then ID if GPL then dex
            % reference
            if isfield(obj, {'calls'}) && ~isempty(obj.calls)
                disp('GPL Calls Detected, array ID variable being filled by dex')
                IDval = obj.arrivalTable.dex+1; % Bug in GPL code leave this in until recieved fix from Tyler XXX
            else
                IDval = obj.arrivalTable.ID;
            end
            
            
            array = [obj.arrivalTable.ArrivalSec(:,...
                [obj.array_struct.master, obj.array_struct.slave(obj.child_idx)])...
                obj.arrivalTable.Location ...
                IDval];
            
            
            % Remove calls that weren't detected on the parent hydrophone
            array = array(~isnan(array(:,1)),:);
            
            % At least 1 tdoa value needed, therefore remove any calls
            % (master) where there aren't at least two arrivals
            array = array(sum(~isnan(array(:,1:3)),2)>= 2,:);
            
            
            
            
            % Sort by start time on the parent hydrophone
            [~,sortidx] = sort(array(:,1));
            array = array(sortidx,:);
            
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Create the Associations and TDOA values %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            % If we are adding in random missassociation then modify the
            % arrival array
            if obj.randomMiss == 1
                
                % Create a matrix that incorporates miss association
                ArrivalsMiss = zeros(size(array));
                ArrivalsMiss(:,1)  = array(:,1);
                ArrivalsMiss(:,end)  = array(:,end);
                
                % step through each child hydrophone and caluculate the
                % TDOA and maximum expected TDOA based on the hydrophone
                % array
                for channum =1:length(obj.child_idx)
                    
                    % Shift the arrivals to simulate clock drift by either
                    % positive or negative values
                    clock_shift = (mod(channum, 2)*2-1) * obj.drift;
                    
                    array(:,channum+1) = array(:,channum+1) + clock_shift;
                    
                    % Allow for Mis Association %%%%%%%%%%%%%%%%%
                    
                    % Distance between two sensors - maximum exptected TDOA
                    % depth between calling whale and array
                    depth_range = obj.hydrophone_struct(obj.array_struct.master).depth...
                        - obj.hydrophone_struct(obj.array_struct.slave(channum)).depth ;
                    
                    % Calculate the horizontal distance between the calling
                    % locations and the hydrophone
                    horizontal_distance = arrayfun(@(lats, lons)...
                        vdist(lats, lons,...
                        obj.hydrophone_struct(obj.array_struct.master).location(1),...
                        obj.hydrophone_struct(obj.array_struct.master).location(2)),...
                        obj.hydrophone_struct(obj.array_struct.slave(channum)).location(1),...
                        obj.hydrophone_struct(obj.array_struct.slave(channum)).location(2));
                    
                    
                    % Maximum expected tdoa between parent and child
                    % (channum)
                    MaxTOA_1 = sqrt(depth_range^2 + ...
                        horizontal_distance.^2)/obj.c;
                    
                    
                    % Step through the arrivals and select a call at random
                    % from all calls within the MaxTOA_1
                    for ii = 1:length(ArrivalsMiss)
                        
                        % Get index of all calls that fall within the
                        % association zone plus the wiggle room
                        corr_idxs = find(...
                            abs(array(:,(channum+1))-array(ii,1))<...
                            MaxTOA_1+obj.assSec);
                        
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
            
            obj.arrivalArray =array;
        end
        
        %% Create table of arrival times (step 1)
        function UpdateArrTable(obj)
            % Use the TDOA table to create a matrix with TDOAs that are
            % clipped based on distance to receiver
            
            % Make life a little easier
            spaceWhaleloc = obj.spaceWhale;
            
            % Create the arrivals table (at) from the first agent
            at = table(...
                spaceWhaleloc.agent(1).call_times', ...
                spaceWhaleloc.agent(1).Arrival_times,...
                spaceWhaleloc.agent(1).RangeKm,...
                spaceWhaleloc.agent(1).location(spaceWhaleloc.agent(1).call_times,:),...
                ones(length(spaceWhaleloc.agent(1).call_times),1),...
                struct2array( spaceWhaleloc.agent(1).tdoa(obj.array_struct.master)),...
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
                    struct2array( spaceWhaleloc.agent(ii).tdoa(obj.array_struct.master)),...
                    'VariableNames',{'callTimes','ArrivalSec', 'RangeKm',...
                    'Location','ID', 'TDOA'});
                at = [at ; Tnew];
            end
            
            % annotate out any detections greater than Xkm from the receiver
            distance_bool = ones(size(at.RangeKm));
            distance_bool(find(at.RangeKm>obj.truncateKm)) = nan;
            
            at.RangeKm = at.RangeKm.*distance_bool;
            at.ArrivalSec = at.ArrivalSec.*distance_bool;
            obj.arrivalTable =at;
        end
        
        %% Clear calculated values (simmat, projections, TODA values, etc)
        function clearCalcValues(obj)
            
            % Clear out all values that need to be calculated leaving the
            % parameters intact. Useful for running multiple experiments
            % using the same agent-setup but different correlation or gap
            % parameters
            obj.arrivalTable=[];
            obj.TDOA_vals=[];

            obj.Cluster_id=[];
            obj.AdjRand=[];
            obj.Sim_mat =[];
            obj.arrivalArray =[];
            obj.chains =[];
            obj.AdjRand=[];
            
            
        end
        
        %% Function for calculating the averageed (across hydrophones)
        % probloc for a given call- NOT normalized or drawn from normal
        % distribution
        function averageLklhd_space = getAvLkHdSpace(obj, callIdx)
            % Returns an averaged likelihood space for a given TDOA
            
            % Pick the set of call delays (number of calls)
            delays = (obj.TDOA_vals(callIdx, :));
            
            
            child_idx = obj.child_idx;
            
            child_hyds = obj.array_struct.slave(obj.child_idx);
            
            
            averageLklhd_space = [];
            
            for ii=1:length(child_idx)
                
                % Get the expected delay spaces for the hydrophone pairs
                toa_space = cell2mat(obj.array_struct.toa_diff(child_idx(ii)+1));
                
                % Delta TOA space
                averageLklhd_space(:,:,ii) = (toa_space - delays(ii));
                
%                 
%                 figure;
%                 contourf(obj.array_struct.longrid, obj.array_struct.latgrid, (averageLklhd_space(:,:,ii)))
%                 colorbar;
%                 caxis([-3 3])
%                 hold on
%                 scatter(obj.arrivalArray(callIdx,end-1), obj.arrivalArray(callIdx,end-2))
            end
            
        end
        
        %% Retruns Normalized Likelihood
        function averageLklhd_space = getTruHdSpace(obj, idx, sig_tot)
            
            % Get the observed TDOA space and normalize
            averageLklhd_space = getAvLkHdSpace(obj, idx);
            
            % set up the vectorization for the PDF
            sigma = ones(size(averageLklhd_space)).*sig_tot;
            
            % Create ambiguity surface and normalize
            averageLklhd_space = normpdf(averageLklhd_space, 0, sigma)./...
                normpdf(0, 0, sigma);
            
            
            if ndims(averageLklhd_space)>1
                
                % sum along third axis, will be normalized later
                averageLklhd_space = nanmean(averageLklhd_space,3);
                
            end
            
            
            
        end
        
        %% Function for reporting the adjusted Rand indices
        function getRand(obj)
            
            % Check if cluster id's are present if not update it
            if isempty(obj.Cluster_id)
                updateClusterID(obj)
            end
            
            
            % Pull out the arrival array for easy reading
            arrival_array = obj.arrivalArray;
            trueClusters =unique(arrival_array(:,end));
            
            % Create smaller clusters based on the cut time
            timeClusters = tempClust(obj);
            
            
            % Create a dummy varaible for the clusters
            NewTruth = (arrival_array(:,end));
            
            
            % step through each true cluster and determine if it crosses a
            % break. If so it gets a new name;
            if obj.tweakedRand ==1
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
            else
                NewTruth = (arrival_array(:,end));
            end
            
            
            
            
            
            % If the arrival array isn't empty, grab the adjusted rand idx
            if ~isempty(arrival_array)
                
                % Allign the true clusters and the predicted clusters
                newClusterId = allignclusters(NewTruth,...
                    obj.Cluster_id)+1;
                
                % Get the adjusted rand (third party)
                [f,~,~,~] = RandIndex(newClusterId, NewTruth);
                obj.AdjRand = f;
            else
                % Otherwise, there were too few data to cluster
                obj.AdjRand = nan;
                disp('No rand values available')
            end
            
            
        end
        
        
        %% Create chains based on time of arrival only
        function TOAChains(obj)
            
            % Check if cluster id's are present if not update it
            if isempty(obj.Cluster_id)
                updateClusterID(obj)
            end
            % Pull out the arrival array for easy reading
            arrival_array = obj.arrivalArray;
            
            % If the arrival array isn't empty, grab the adjusted rand idx
            if ~isempty(arrival_array)
                
                % Allign the true clusters and the predicted clusters
                newClusterId = allignclusters(arrival_array(:,end),...
                    obj.Cluster_id);
                
                % Get the adjusted rand (third party)
                [~,f,~,~] = RandIndex(newClusterId, arrival_array(:,end));
                obj.AdjRand = f;
            else
                % Otherwise, there were too few data to cluster
                obj.AdjRand = nan;
            end
            
            
        end
        
        %% Function for creating predicted cluster ID's for the agents
        function updateClusterID(obj)
            
            % If chains have not been updated, do so
            if isempty(obj.chains)
                disp('Updating Simulation Matrix')
                updateChains(obj)
            end
            
            % Simply create an array with the predicted clusterid
            Cluster_id = zeros(length(obj.arrivalArray),1)+length(obj.chains)+1;
            
            % Step through each chain and grab the cluster
            for ii=1:length(obj.chains);
                Cluster_id(obj.chains(ii).index) = ii;
                
            end
            
            obj.Cluster_id =Cluster_id;
        end
        
        %% Function for creating the cluster chains (aka ladder linkages)
        function updateChains(obj)
            
            % Original design via Glen Lerley
            if isempty(obj.Sim_mat)
                disp('Updating Simulation Matrix')
                UpdateSimMat(obj)
            end
            
            
            % Ladder linkages for new reduced plot
            cutoff =obj.cutoff; % Correlation threshold cut off
            time_cut = obj.time_cut;
            
            time = obj.arrivalArray(:,1);
            
            % Elapsed time since first call
            time=(time-time(1));  % seconds
            time0=time; % keep original copy for later indexing
            
            % Replicate Corr_coef_map for shrinking matrix
            Corr_coef_map = obj.Sim_mat;
            [s1,s2]=size(Corr_coef_map);
            nc=1; % cluster counter
            clear chain dex ss
            
            % While there are still rows in the matrix
            while s1>1
                flag=1;
                
                % call element index, start at 1 and examine columnsmaxas
                elt=1; % first element of cluster is always upper left
                % hand element
                index=[];
                cluster=[];
                
                
                % find clusters exceeding the threshold and not exceeding
                % the time cutoff
                
                while flag == 1
                    % Add call to current cluster, retaining the time and index
                    % of the call element
                    cluster = [cluster, time(elt)];
                    index   = [index, elt];
                    
                    % entire column for each element
                    vec=Corr_coef_map(:,elt);
                    
                    % Only look at the elements below
                    vec=vec(elt+1:end);
                    
                    % Look for correlations that exceed the threshold
                    k=find(vec>cutoff);
                    
                    if isempty(k)
                        % Nobody met the criterion
                        flag=0;
                    else
                        % Set next call to the first that meets criterion
                        elt=k(1)+elt;
                    end
                end % flag test
                
                % now cut cluster at first large jump:
                k=find(diff(cluster)>time_cut);
                
                % If any objects in the cluster
                if length(k)>0
                    k=k(1);
                    % Remove everything after the first large gap
                    cluster=cluster(1:k);
                    index=index(1:k);
                    
                    % Remove calls with separation times of less than one
                    % second
                    if length(cluster)>1
                        cluster = [cluster(1) cluster(find(diff(cluster)>1))];
                        index=[index(2) index(find(diff(cluster)>1))];
                        cluster = cluster(cluster>0);
                        index = index(index>0);
                    end
                end
                
                % aggregate clusters
                chain(nc).links=cluster;
                chain(nc).n = length(cluster);
                for j=1:length(cluster)
                    chain(nc).index(j)=find(time0==cluster(j));
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
        
        %% Compare two probability grid spaces (Help me!)
        function simValue = compareLklhdSpace(obj,...
                ProjectedLklhdSpace, nextLklhdSpace)
            
            
            
            % more of the same
            
            
            nextArea = find(nextLklhdSpace>0.01);
            projArea = find(ProjectedLklhdSpace>.01);
            
            commonvals = intersect(nextArea, projArea);
            
            simValue = sum(nextLklhdSpace(commonvals))/sum(sum(nextLklhdSpace));



            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
           
        %% Draw the Agents
        function drawAgents(obj)
            
            % Update the arrival array
            if isempty(obj.arrivalArray)
                UpdateArrArray(obj);
            end
            
            % Hydrophone locations
            hyd_table = struct2table(obj.hydrophone_struct);
            
            TimeColorVals = parula(obj.spaceWhale.param_sim.dur+2);
            ColorVals = lines( max([length(obj.spaceWhale.agent), max(obj.Cluster_id)]));
            child_hyds = obj.array_struct.slave(obj.child_idx);
            
            
            figure;
            
            subplot(3,1,2)
            hold on
            scatter(obj.arrivalArray(:,end-1),...
                obj.arrivalArray(:,end-2),[],...
                ColorVals(obj.arrivalArray(:,end),:), 'f')
            
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
            
            subplot(3,1,1)
            
            hold on
            scatter(obj.arrivalArray(:,end-1),...
                obj.arrivalArray(:,end-2),[],...
                [TimeColorVals(round(obj.arrivalArray(:,1)),:)],...
                'f')
            
            scatter(hyd_table.location(:,2), hyd_table.location(:,1), 80, 'k', 'filled', 'd')
            scatter(...
                hyd_table.location([obj.array_struct.master, obj.array_struct.slave(obj.child_idx)],2),...
                hyd_table.location([obj.array_struct.master, obj.array_struct.slave(obj.child_idx)],1),...
                'r', 'filled', 'd')
            
            colorbar
            title('True Cluster')
            
            if ~isempty(obj.Cluster_id)
                
                subplot(3,1,3)
                hold on
                scatter(obj.arrivalArray(:,end-1), obj.arrivalArray(:,end-2),[],...
                    ColorVals(obj.Cluster_id,:), 'f')
                
                scatter(hyd_table.location(:,2), hyd_table.location(:,1),...
                    80, 'k', 'filled', 'd')
                
                scatter(...
                    hyd_table.location([obj.array_struct.master, obj.array_struct.slave(obj.child_idx)],2),...
                    hyd_table.location([obj.array_struct.master, obj.array_struct.slave(obj.child_idx)],1),...
                    'r', 'filled', 'd')
                
                titlestr = [obj.titleStr, ' Expected Clusters'];
                title(titlestr)
                
            end
            
            
            
        end
        
        %% Draw the simulation matrix
        function drawSimMat(obj)
            figure;
            [nr,nc] = size(obj.Sim_mat);
            h =pcolor(flipud((obj.Sim_mat)));
            set(h, 'EdgeColor', 'none');
            ax = gca;
            ax.YTickLabel = flipud(ax.YTickLabel)
            colormap(jet);
            colorbar
            title(obj.titleStr)
            xlabel('Call ID')
            ylabel('Call ID')
            
            
        end
    end
    
end
