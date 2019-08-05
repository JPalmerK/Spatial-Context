function [propCorrect, NClusters, nCorrect] = RunComp(examp, parent, hyd, truth, corr_thresh)

% Input -
% Examp: simulation class object pre-processed for similarity matrix
% Truth: Validated Raven selection table
% Threshold: correlation threshold
%% Run the ladder linkages with the given correlation threshold

if ~strcmp(examp.titleStr,   'Baseline - toa Only')
examp.cutoff = corr_thresh;
updateClusterID(examp)
end




localize_struct = examp.localize_struct;
arivalArr = examp.arrivalArray;
%% Create a raven table from the detections
RavenTable = array2table(zeros(0,12),...
    'VariableNames',{'Selection', 'View', 'Channel',...
    'BeginS', 'EndS', 'LowF', 'HighF', 'Lat', 'Lon','ClusterId',...
    'MtlbDtStr','nPairs'});

% Loop through the child indexes and create the arrival tables

hyds = [parent, localize_struct.hyd(parent).array.slave(examp.child_idx)];
% Calls

hyd_idx = [parent  localize_struct.hyd(parent).array.slave];


calls =struct2table(hyd(parent).detection.calls);

calls.duration = (calls.end_time-calls.start_time)/2000;
child_idx =1;

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


gpldat =  RavenTable;
gpldat.mstart = datenum('20090328', 'yyyymmdd')+gpldat.BeginS/86400;
gpldat.mend =datenum('20090328', 'yyyymmdd')+gpldat.EndS/86400;
gpldat.spp = repmat({'unknown'}, [height(gpldat),1]);
gpldata.Lat =RavenTable.Lat;
gpldat.Lon = RavenTable.Lon;

% step through and add the species where available
for ii=1:height(gpldat)
    
    gpl_time = (gpldat.mstart(ii)+gpldat.mend(ii))/2;
    
    
    % Canidate calls
    linkedID = (...
        intersect(...
        find(truth.Channel==gpldat.Channel(ii)),...
        find((abs(truth.mid - gpl_time)*24*60*60)<.25)));
    
    
    if length(linkedID)>=1
        spp =truth.Species(linkedID);
        gpldat.spp(ii) = spp(1);
        %disp(num2str(ii));
    end
    
    
    
    
end



h1 = gpldat(logical(strcmp(gpldat.spp, 'hb')),:);
h2 = gpldat(logical(strcmp(gpldat.spp, 'rw')),:);
h3 = gpldat(logical(~strcmp(gpldat.spp, 'unknown')),:);
%
% figure;
% subplot(2,1,1)
% title('Truth: Known Species IDs')
% hold on;
% scatter(h1.Lon, h1.Lat,[], 'f', 'b')
% scatter(h2.Lon, h2.Lat,[], 'f', 'r')
% %scatter(hyd_arr(:,2),hyd_arr(:,1), 80, 'k', 'filled', 'd')
% legend('Humpback', 'Right Whale')
%
% subplot(2,1,2)
% title('Clustered Results')
% clust_ids = unique(h3.ClusterId);
% colors = jet(length(clust_ids));
% colors=colors(randsample(1:length(colors),length(colors)),:);

% for ii=1:length(clust_ids)
%     hold on
%     ha = h3(h3.ClusterId==clust_ids(ii),:);
%     scatter(ha.Lon, ha.Lat,[], colors(ii,:), 'f');
%     
% end

N = (unique(h3.ClusterId));
% legendCell = char(transpose(cellstr(num2str(N, '%-d'))));
% legend(legendCell);

%scatter(hyd_arr(:,2),hyd_arr(:,1), 80, 'k', 'filled', 'd')





k2 = ~strcmp(gpldat.spp, 'unknown');

wrongclust = h3((h3.ClusterId==1),:);


%gpldatout = gpldat(gpldat.Channel==5,:);
gpldatout = gpldat(k2,:);
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

NClusters = length(unique(examp.Cluster_id));



end











