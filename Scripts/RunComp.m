function [propCorrect, NClusters, nCorrect,GPLSummary] = RunComp(examp, parent, hyd, truth, corr_thresh)

% Input -
% Examp: simulation class object pre-processed for similarity matrix
% Truth: Validated Raven selection table
% Threshold: correlation threshold





localize_struct = examp.localize_struct;
arivalArr = examp.arrivalArray;
%% Create a raven table from the detections
RavenTable = array2table(zeros(0,12),...
    'VariableNames',{'Selection', 'View', 'Channel',...
    'BeginS', 'EndS', 'LowF', 'HighF', 'Lat', 'Lon','ClusterId',...
    'MtlbDtStr','nPairs'});

% Loop through the child indexes and create the arrival tables
hyd_idx = [parent  localize_struct.hyd(parent).array.slave];


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
    end_times = start_times + calls.duration(localize_struct.hyd(parent).dex);
    
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
    
    %     % Index of the parent call
    %     call_ids =examp.localize_struct.hyd(parent).dex(arra_ids)-1;
    
    
    
    
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
        ClusterId,...
        matlabDate,...
        narrivals,...
        'VariableNames',{'Selection', 'View', 'Channel',...
        'BeginS', 'EndS', 'LowF', 'HighF', 'Lat', 'Lon','ClusterId',...
        'MtlbDtStr','nPairs'});
    
    RavenTable =[RavenTable; aa];
    
end

idxTruth = find(truth.LowFreq_Hz_ < hyd(5).detection.parm.freq_hi);
truth_trimmed = truth(idxTruth,:);

RavenTable.mstart = datenum('20090328', 'yyyymmdd')+RavenTable.BeginS/86400;
RavenTable.mend =datenum('20090328', 'yyyymmdd')+RavenTable.EndS/86400;
RavenTable.spp = repmat({'unknown'}, [height(RavenTable),1]);
RavenTable.Lat =RavenTable.Lat;
RavenTable.Lon = RavenTable.Lon;

wiggleroom_days = .5/60/60/24;
% step through and add the species where available
for ii=1:height(RavenTable)
    
    gpl_time = (RavenTable.mstart(ii)+RavenTable.mend(ii))/2;
    bb = truth_trimmed.Channel==RavenTable.Channel(ii);
    cc = truth_trimmed.mend> (gpl_time-wiggleroom_days);
    dd = truth_trimmed.mstart< (gpl_time +wiggleroom_days);
    
    linkedID = find(prod([bb, cc, dd],2));
    
    if length(linkedID)>=1
        spp =truth_trimmed.Species(linkedID);
        RavenTable.spp(ii) = spp(1);
        %disp(num2str(ii));
    end
    
    clear bb cc dd
    
    
end



k2 = ~strcmp(RavenTable.spp, 'unknown');


gpldatout = RavenTable(k2,:);
gpldatout.voting=repmat({'Unk'},[height(gpldatout),1]);
cluster_ids = unique(gpldatout.ClusterId);

for ii = 1:length(cluster_ids)
    
    clus = cluster_ids(ii);
    spp = gpldatout.spp(gpldatout.ClusterId==clus);
    
    spp = majorityvote(spp);
    gpldatout.voting(gpldatout.ClusterId==clus) = spp;
    
    
end

% Number of right whales and number of humpabcks
rw_sum = sum(cellfun(@strcmp, gpldatout.spp, repmat({'rw'},height(gpldatout),1)));
hb_sum = sum(cellfun(@strcmp, gpldatout.spp, repmat({'hb'},height(gpldatout),1)));


% Proportion of correct answers
correct = cellfun(@strcmp, gpldatout.spp, gpldatout.voting);
nCorrect = sum(correct);
propCorrect = nCorrect/height(gpldatout);

NClusters = length((gpldatout.ClusterId))/length(unique(gpldatout.ClusterId));





end







