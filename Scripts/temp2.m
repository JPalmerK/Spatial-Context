
calls =hyd(5).detection.calls

% Find GPL indexes of the calls that coincide with the truth tables

close all; 
    % Populate data and parameters
    examp = clusterClassGPL();
    examp.array_struct = array_struct;
    examp.hydrophone_struct = hydrophone_struct;
    examp.cutoff = .75;
    examp.time_cut = 70*60;
    examp.randomMiss =0;
    examp.s =8;
    examp.child_idx = [1:8];
    examp.localize_struct =localize_struct;
    examp.limitTime =24*60*60;
    examp.maxEltTime =60*10;
    examp.calls =localize_struct.hyd(5).calls;
    examp.callParm =  hyd(5).detection.parm;       
    
    examp.clearCalcValues
    UpdateArrTable(examp) % Run this first!
    UpdateArrArray(examp)
    
parent = array_struct.master;
arivalArr = examp.arrivalArray;

RavenTable = array2table(zeros(0,9),...
    'VariableNames',{'Selection', 'View', 'Channel',...
    'BeginS', 'EndS', 'LowF', 'HighF','ClusterId', 'MtlbDtStr'});

% Loop through the child indexes and create the arrival tables
hyds = [array_struct.master, array_struct.slave(examp.child_idx)];

% Calls
calls = struct2table(calls);

% frequency information
f_low =  hyd(5).detection.parm.freq_lo; % frequency limits
f_high = hyd(5).detection.parm.freq_hi;

df = (f_high-f_low)/(examp.callParm.bin_hi-examp.callParm.bin_lo);
for ii =1:length(hyds)
    n_calls =sum(~isnan(arivalArr(:,ii)));
    call_ids =examp.localize_struct.hyd(parent).dex(logical(~isnan(arivalArr(:,ii))));
    
    % Error in 
    call_ids = call_ids;
    
    % Start Time
    start_times = calls.start_time(call_ids)/examp.fs;
    
    % End Times (DOOoooOOOoOOOoOOM!)
    end_times = calls.end_time(call_ids)/examp.fs;
    
    % Get the high and low frewuency based on the spectrogram
    % parameters
    low_f = zeros(length(end_times),1);
    high_f = low_f;
    
    for jj=1:n_calls
        [rr, cc]=ind2sub(calls.cm(call_ids(jj)).size, calls.cm(call_ids(jj)).index);
        low_f(jj) = f_low+(df*min(rr));
        high_f(jj) = f_low+(df*max(rr));
    end
    
    
    Selection =  [height(RavenTable)+1:height(RavenTable)+n_calls]';
    View =repmat({'Spectrogram'},[n_calls,1]);
    Channel =  repmat(hyds(ii), [n_calls,1]);
    ClusterId = ones([n_calls, 1]);
    matlabDate = calls.julian_start_time(logical(~isnan(arivalArr(:,ii))));
    
    aa = table(Selection,...
        View, ...
        Channel,...
        start_times,...
        end_times,...
        low_f,...
        high_f,...
        ClusterId,...
        matlabDate,...
        'VariableNames',{'Selection', 'View', 'Channel',...
        'BeginS', 'EndS', 'LowF', 'HighF', 'ClusterId','MtlbDtStr'});
    
    RavenTable =[RavenTable; aa];
    
end


% Find the indicies that coinced with the truth table
% Load the calls
truth = readtable('C:\Users\Kaitlin\Desktop\NOPP6_EST_20090328.selections.txt');
truth.mstart = datenum('20090328', 'yyyymmdd')+truth.BeginTime_s_/60/60/24;
truth.mend = datenum('20090328', 'yyyymmdd')+truth.EndTime_s_/60/60/24;

gpldat =  RavenTable;
gpldat.mstart = gpldat.MtlbDtStr;
gpldat.mend =gpldat.mstart+(gpldat.EndS-gpldat.BeginS)/2/60/60/24;
gpldat.spp = repmat({'unknown'}, [height(gpldat),1]);


% step through and add the species where available
for ii=1:height(gpldat)
    
    gpl_time = (gpldat.mstart(ii)+gpldat.mend(ii))/2;
    gpl_f = (gpldat.LowF(ii)+gpldat.HighF(ii))/2;

    % Canidate calls
    linkedID = find(truth.Channel==gpldat.Channel(ii) &...
        truth.mstart<= gpl_time &...
        truth.mend>= gpl_time &...
        truth.LowFreq_Hz_<= gpl_f &...
        truth.HighFreq_Hz_>= gpl_f);
    if length(linkedID)==1
    spp =truth.Species(linkedID)
    gpldat.spp(ii) = spp;
    
    elseif length(linkedID)>1
       print('blarg')
    end

    

end

k2 = find(~strcmp(gpldat.spp, 'unknown'));
clear gpldat











