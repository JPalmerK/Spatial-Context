function GPLSummary = GPLRavenSummaryComp(GPLdat, hyd, truth, thresh, parent)

% Compare the number of detections annotated by the analyst and the number
% of detections that match from GPL. GPL can either be the raw number of
% detections (hyd) or the output from the localizeation algorithim
% (localize_struct)

% Determine if GPL input is hydrophone struct or localize struct
if size(GPLdat,2)>1
    % hydrophone struct
    GPLdat = trimCallsCorrThresh(GPLdat, thresh);
    
    

    flag =0;
    
else
    %localize struct, pull the hydrophone data
    GPLdat = trimlocalize_struct(GPLdat,...
    hyd, thresh);
  
    calls =struct2table(hyd(parent).detection.calls);
    calls.duration = (calls.end_time-calls.start_time)/2000;
    child_idx =1;
    % Loop through the child indexes and create the arrival tables
    hyd_idx = [parent  GPLdat.hyd(parent).array.slave];
    flag =1;
    child_idx =1;

end




RavenTable = table();

% Create the raven table from the localize structure
for ii =1:length(GPLdat.hyd)
    
    
    % Current hydrophone id
    hyd_id = ii;
    
    % Delay time to add, if parent add nothing
    if (hyd_id~= parent) & flag
        hyd_delay = GPLdat.hyd(parent).delays(:,child_idx);
        child_idx= child_idx+1;
    else
        hyd_delay =0;
    end
    
    % Arrival times
    Arrival_times = GPLdat.hyd(parent).rtimes'/2000+hyd_delay; 
    start_times = Arrival_times-(0.5*calls.duration(GPLdat.hyd(parent).dex));
    end_times = start_times + calls.duration(GPLdat.hyd(parent).dex);
    
    % Low and High FRequency
    low_f = calls.flow(GPLdat.hyd(parent).dex);
    high_f =calls.high(GPLdat.hyd(parent).dex);
    
    % Lat and Lon
    lat = squeeze(GPLdat.hyd(parent).coordinates(end,1,:));
    lon = squeeze(GPLdat.hyd(parent).coordinates(end,2,:));
    
    % Call detected on the hydrophone
    arra_ids = find(~isnan(GPLdat.hyd(parent).cross_score(:,ii)));
    n_calls =length(arra_ids);
    narrivals = sum(~isnan(GPLdat.hyd(parent).cross_score),2);
    
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






% number of detections by channel
nGPLdet = zeros([1,9]);
nGPLMN = nGPLdet;
nGPLEG = nGPLdet;
nGPLFP = nGPLdet;
nMNTP = nGPLdet;
nEGTP = nGPLdet;

for ii =1:9
    %Number of GPL detections, total per channel
    nGPLdet(ii) = sum(RavenTable.Channel == chan_ids(ii));
    
    
    sppBInary = cellfun(@strcmp, RavenTable.spp, repmat({'hb'},height(RavenTable),1));
    ChannBinary = RavenTable.Channel == chan_ids(ii);
    % GPL that were associated with humpback per channel
    nGPLMN(ii) = sum(and(sppBInary, ChannBinary));
    
    sppBInary = cellfun(@strcmp, RavenTable.spp, repmat({'rw'},height(RavenTable),1));
    ChannBinary = RavenTable.Channel == chan_ids(ii);
    % GPL that were associated with rightwhale per channel
    nGPLEG(ii) = sum(and(sppBInary, ChannBinary));


    sppBInary = cellfun(@strcmp, RavenTable.spp, repmat({'unknown'},height(RavenTable),1));
    ChannBinary = RavenTable.Channel == chan_ids(ii);
    %Total GPL that were false positve
    nGPLFP(ii) = sum(and(sppBInary, ChannBinary));
    
    
    sppBInary = cellfun(@strcmp, truth_trimmed.Species, repmat({'hb'},height(truth_trimmed),1));
    ChannBinary = truth_trimmed.Channel == chan_ids(ii);
    nMNTP(ii) = sum(and(sppBInary, ChannBinary));
    
     sppBInary = cellfun(@strcmp, truth_trimmed.Species, repmat({'rw'},height(truth_trimmed),1));
    ChannBinary = truth_trimmed.Channel == chan_ids(ii);
    nEGTP(ii) = sum(and(sppBInary, ChannBinary));

end



GPLSummary = table(nGPLdet', nGPLMN', nMNTP', nGPLEG', nEGTP', nGPLFP',...
    'Variable',{'GPLDet', 'GPLDetMn', 'TPMN','GPLDetEg','TPEG', 'nGPLFP'});
GPLSummary.Channel = chan_ids;

