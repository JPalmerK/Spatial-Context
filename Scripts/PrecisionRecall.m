function [Clustered, Detector, RavenTable, FNBaseline] = ...
    PrecisionRecall(examp, parent, hyd, truth)

% Input -
% Examp: simulation class object pre-processed for similarity matrix
% Truth: Validated Raven selection table
% Threshold: correlation threshold



localize_struct = examp.localize_struct;
arivalArr = examp.arrivalArray;
%% Create a raven table from the detections
RavenTable = array2table(zeros(0,13),...
    'VariableNames',{'Selection', 'View', 'Channel',...
    'BeginS', 'EndS', 'LowF', 'HighF', 'Lat', 'Lon','Scores','ClusterId',...
    'MtlbDtStr','nPairs'});

% Loop through the child indexes and create the arrival tables
hyd_idx = [parent  localize_struct.hyd(parent).array.slave];

% Convert localize struct to data table/raven compatable file

    calls =struct2table(hyd(parent).detection.calls);
    calls.duration = (calls.end_time-calls.start_time)/2000;
    child_idx =1;
    
    % Create the raven table from the localize structure
    for ii =1:length(hyd_idx)
        
        
        % Current hydrophone id
        hyd_id = hyd_idx(ii);
        
        % Delay time to add, if parent add nothing
        if hyd_id~= parent
            hyd_delay = localize_struct.hyd(parent).delays(:,child_idx);
            child_idx= child_idx+1;
        else
            hyd_delay =0;
        end
        
        % Arrival times
        Arrival_times = localize_struct.hyd(parent).rtimes'/2000+hyd_delay;
        start_times = Arrival_times-(0.5*calls.duration(localize_struct.hyd(parent).dex));
        end_times = start_times + (0.5*calls.duration(localize_struct.hyd(parent).dex));
        
        % Low and High FRequency
        low_f = calls.flow(localize_struct.hyd(parent).dex);
        high_f =calls.high(localize_struct.hyd(parent).dex);
        
        % Lat and Lon
        lat = squeeze(localize_struct.hyd(parent).coordinates(end,1,:));
        lon = squeeze(localize_struct.hyd(parent).coordinates(end,2,:));
        
        % Call detected on the hydrophone
        arra_ids = find(~isnan(arivalArr(:,ii)));
        n_calls =length(arra_ids);
        narrivals = sum(~isnan(arivalArr(arra_ids,2:9)),2);
        
        % Detector Scores
        if hyd_id==parent
            scores = localize_struct.hyd(parent).detectorScore(arra_ids);
            scores = reshape(scores,[length(scores),1 ]);
        else
            scores = zeros(size(narrivals))/0;
        end
        
        Selection =  [height(RavenTable)+1:height(RavenTable)+n_calls]';
        View =repmat({'Spectrogram'},[n_calls,1]);
        Channel =  repmat(hyd_idx(ii), [n_calls,1]);
        ClusterId = examp.Cluster_id(arra_ids);
        matlabDate = calls.julian_start_time(arra_ids);
        
        aa = table(Selection,...
            View, ...
            Channel,...
            start_times(arra_ids),...
            end_times(arra_ids),...
            low_f(arra_ids),...
            high_f(arra_ids),...
            lat(arra_ids),...
            lon(arra_ids),...
            scores,...
            ClusterId,...
            matlabDate,...
            narrivals,...
            'VariableNames',{'Selection', 'View', 'Channel',...
            'BeginS', 'EndS', 'LowF', 'HighF', 'Lat',...
            'Lon','Scores','ClusterId',...
            'MtlbDtStr','nPairs'});
        
        RavenTable =[RavenTable; aa];
        
    end
    
    % Trim the turth table by the frequency limits of the detector (don't
    % penelize calls outside of the detector frequency range)
    idxTruth = find(truth.LowFreq_Hz_< hyd(parent).detection.parm.freq_hi);
    truth_trimmed = truth(idxTruth,:);
    truth_trimmed = truth_trimmed(truth_trimmed.Channel==parent,:);
    truth_trimmed.Sec = (truth_trimmed.mstart-min(truth_trimmed.mstart))*60*60*24;
    truth_temp = truth_trimmed;
    
    % Add some variables to the GPL Raven table for comparison
    RavenTable = RavenTable(RavenTable.Channel == parent,:);
    RavenTable.mstart = datenum('20090328', 'yyyymmdd')+...
        RavenTable.BeginS/86400;
    RavenTable.mend =datenum('20090328', 'yyyymmdd')+RavenTable.EndS/86400;
    RavenTable.spp = repmat({'unknown'}, [height(RavenTable),1]);
    
    
    %Trim to periods where calls interleave
    RavenTable.Sec = (RavenTable.mstart-floor(RavenTable.mstart))*60*60*24;
    RavenTable = RavenTable(RavenTable.Sec>=min(truth_trimmed.BeginTime_s_) &...
        RavenTable.Sec<= max(truth_trimmed.BeginTime_s_),:);
   
    
    % Overlap between the GPL detections and truth
    wiggleroom_days = .3/60/60/24;
    
    % Link truth species classifications with GPL detections. Add species to
    % GPL and knock out truth row such that the FN rate can be accurately
    % accounted for.
    for ii=1:height(RavenTable)
        
        if height(truth_temp)>0
            % Logicals for finding calls overlapping with GPL detections
            gpl_time = (RavenTable.mstart(ii)+RavenTable.mend(ii))/2;
            
            linkedID = find(truth_temp.Channel==RavenTable.Channel(ii) &...
                truth_temp.mend> (gpl_time-wiggleroom_days) &...
                truth_temp.mstart< (gpl_time +wiggleroom_days));
            
            
            if length(linkedID)>=1
                spp =truth_temp.Species(linkedID);
                RavenTable.spp(ii) = spp(1);
                truth_temp(linkedID,:)=[];
            end
        end
        
        
        
    end
    
    RavenTable.voting= zeros([height(RavenTable),1])/0;
    cluster_ids = unique(RavenTable.ClusterId);
    
    for jj = 1:length(cluster_ids)
        
        clus = cluster_ids(jj);
        
        corrScores = RavenTable.Scores(RavenTable.ClusterId==clus);
        corrScores = corrScores(~isnan(corrScores));
        %LR = geomean(corrScores);
        %LR = max(corrScores);
        LR = log(prod(corrScores./(1.00001-corrScores)))/length(corrScores);
        RavenTable.voting(RavenTable.ClusterId==clus)=LR;
        
        
        
    end
    
    % threhsold scores for calculating precision recall
    RavenTable.LRscores = log(RavenTable.Scores./...
        (1.0000001 - RavenTable.Scores));
    
    scores = unique(RavenTable.LRscores);
    
    
    %Pre allocate output
    RecallDetOnly = zeros(1, length(scores)-1);
    PrecisionDetOnly = RecallDetOnly;
    RecallClustered= RecallDetOnly;
    PrecisionClustered = RecallDetOnly;
    ErrorRateClustered = RecallDetOnly;
    ErrorRateDetector = RecallDetOnly;
    
    
    
    % False positives that were missed
    FNBaseline = length(find(strcmp(truth_temp.Species, 'rw')));
    FNBaseline =0;
    
    for jj = 1:length(scores)
        
        % Detector only Precision Recall
        TP(jj) = sum((RavenTable.LRscores >= scores(jj)).*...
                        strcmp(RavenTable.spp, 'rw'));
       
        
        % Labeled right whale but were not right whale
        FP(jj) = sum((RavenTable.LRscores >= scores(jj)).*...
                    ~strcmp(RavenTable.spp, 'rw'));
        
        % Not labeled right whale but were right whale
        FN(jj) = FNBaseline + ...
                    sum((RavenTable.LRscores < scores(jj)).*...
                        strcmp(RavenTable.spp, 'rw'));
        

        
        % Number of detections that were labeled RW and were in agreement with the
        % truth value
        TPv(jj) = sum((RavenTable.voting >= (scores(jj))).*...
                        strcmp(RavenTable.spp, 'rw'));
        
        % Labeled right whale but were not right whale
        FPv(jj) = sum((RavenTable.voting >= (scores(jj))).*...
                        ~strcmp(RavenTable.spp, 'rw'));
        
        % Not labeled right whale but were right whale
        FNv(jj) = FNBaseline+ ...
                    sum((RavenTable.voting <(scores(jj))).*....
                    strcmp(RavenTable.spp, 'rw')); 
        
    end
    
    
    RecallDetOnly = TP./(TP+FN);
    PrecisionDetOnly = TP./(TP+FP);

    
    
    % Add detections that were trimmed via thresholding to the fn's
    RecallClustered= TPv./(TPv+FNv);
    PrecisionClustered = TPv./(TPv+FPv);

    
    Clustered= struct();
    Clustered.Precision = PrecisionClustered;
    Clustered.Recall =RecallClustered;
    
    
    Detector= struct();
    Detector.Precision = PrecisionDetOnly;
    Detector.Recall =RecallDetOnly;
    
    
    
    








