
% Run experiments
close all; clear all
clear classes; clc

whereAmI = loadSimspace();
dclde_2013_meta = xlsread(whereAmI{1});
load(whereAmI{2})
load(whereAmI{3})
load(whereAmI{4})
load(whereAmI{5})
clear whereAmI


hydrophone_struct= struct();
for ii=1:size(dclde_2013_meta,1)
    hydrophone_struct(ii).name = num2str(dclde_2013_meta(ii,1));
    hydrophone_struct(ii).location = dclde_2013_meta(ii,[11:12]);
    hydrophone_struct(ii).depth= abs(dclde_2013_meta(ii, 13));
    hydrophone_struct(ii).channel=ii;
end
hyd_arr = struct2cell(hydrophone_struct);
hyd_arr =vertcat(hyd_arr{2,:,:});


parent =5;
fs = 2000;
ssp = localize_struct.parm.ssp;
grid_depth = localize_struct.parm.grid_depth;


calls_arrivals =struct2table(hyd(5).detection.calls);

% Find GPL indexes of the calls that coincide with the truth tables

close all; 
    % Populate data and parameters
    examp = clusterClassGPL();
    examp.array_struct = array_struct;
    examp.hydrophone_struct = hydrophone_struct;
    examp.cutoff = .75;
    examp.time_cut = 10*60;
    examp.randomMiss =0;
    examp.s =8;
    examp.child_idx = [1:8];
    examp.localize_struct =localize_struct;
    examp.limitTime =3*60*60;
    examp.maxEltTime =60*10;
    examp.calls =calls_arrivals;
    examp.callParm =  hyd(5).detection.parm;       
    
    examp.clearCalcValues
    UpdateArrTable(examp) % Run this first!
    UpdateArrArray(examp)
    examp.simMatTDOAonly()
    %examp.simMatIdealNewSim();
    updateClusterID(examp)
    examp.updateClusterID;
    examp.drawAgents
    
parent = array_struct.master;
arivalArr = examp.arrivalArray;
%% Create a raven table from the detections
RavenTable = array2table(zeros(0,9),...
    'VariableNames',{'Selection', 'View', 'Channel',...
    'BeginS', 'EndS', 'LowF', 'HighF','ClusterId', 'MtlbDtStr'});

% Loop through the child indexes and create the arrival tables
hyds = [array_struct.master, array_struct.slave(examp.child_idx)];

% Calls

hyd_idx = [array_struct.master array_struct.slave]


calls =struct2table(hyd(5).detection.calls);

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
    
    start_times = Arrival_times-(0.5*calls.duration(localize_struct.hyd(parent).dex+1));
    end_times = start_times + calls.duration(localize_struct.hyd(parent).dex+1);
    
    % Low and High FRequency
    low_f = calls.flow(localize_struct.hyd(parent).dex+1);    
    high_f =calls.high(localize_struct.hyd(parent).dex+1); 
    
    
    % Call detected on the hydrophone
    arra_ids = find(~isnan(arivalArr(:,ii)));
    n_calls =length(arra_ids);
    
    % Index of the parent call
    call_ids =examp.localize_struct.hyd(parent).dex(arra_ids)-1;    
    
    
    
    
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
        ClusterId,...
        matlabDate,...
        'VariableNames',{'Selection', 'View', 'Channel',...
        'BeginS', 'EndS', 'LowF', 'HighF', 'ClusterId','MtlbDtStr'});
    
    RavenTable =[RavenTable; aa];
    
end
% 
% fname ='/home/kpalmer/Desktop/cluster_test.txt';
% writetable(RavenTable(:,[1:9, end]), fname, 'Delimiter', '\t', 'WriteVariableNames', false)    
% header=strcat('Selection\tView\tChannel\tBegin Time (s)\tEnd Time (s)',...
%     '\tLow Freq (Hz)\tHigh Freq (Hz)\',...
%     'tClusterID\tMT1\tSpp\n');
% 
% 
% % export the file
% % - Read input file data content (without header).
%  fid_in = fopen(fname, 'r') ;             % Open input file for reading.
%  fgetl(fid_in) ;                                 % Skip header line in input file.
%  content = fread(fid_in) ;                       % Read rest of input file.
%  fclose(fid_in) ;                                % Close input file.
%  
%  
%  % - Write output file: new header + previous data content.
%  fid_out = fopen(fname, 'w') ;  % Open input file for writing.
%  fprintf(fid_out, header') ;       % Output new header.
%  fwrite(fid_out, content) ;                      % Output previous data content.
%  fclose(fid_out) ;  
% 
% 
% 
% disp('done')

%% Compare GPL raven table to truth raven table
% Find the indicies that coinced with the truth table
% Load the calls
%truth = readtable('C:\Users\Kaitlin\Desktop\NOPP6_EST_20090328.selections.txt');
truth = readtable('/home/kpalmer/Desktop/NOPP6_EST_20090328.selections.txt');
truth.mstart = datenum('20090328', 'yyyymmdd')+truth.BeginTime_s_/86400;
truth.mend = datenum('20090328', 'yyyymmdd')+truth.EndTime_s_/86400;
truth.mid = (truth.mend+truth.mstart)/2;

gpldat =  RavenTable;
gpldat.mstart = datenum('20090328', 'yyyymmdd')+gpldat.BeginS/86400;
gpldat.mend =datenum('20090328', 'yyyymmdd')+gpldat.EndS/86400;
gpldat.spp = repmat({'unknown'}, [height(gpldat),1]);


% step through and add the species where available
for ii=1:height(gpldat)
    
    gpl_time = (gpldat.mstart(ii)+gpldat.mend(ii))/2;


    % Canidate calls
    linkedID = (...
        intersect(...
        find(truth.Channel==gpldat.Channel(ii)),...
        find((abs(truth.mid - gpl_time)*24*60*60)<.25)));
    
    
    if length(linkedID)>=1
    spp =truth.Species(linkedID)
    gpldat.spp(ii) = spp;
    disp(num2str(ii))
    end
    

    

end

k2 = find(~strcmp(gpldat.spp, 'unknown'));




%%
fname = '/home/kpalmer/Desktop/clustercomp.txt'

header=strcat('Selection\tView\tChannel\tBegin Time (s)\tEnd Time (s)',...
    '\tLow Freq (Hz)\tHigh Freq (Hz)\',...
    'tClusterID\tMT1\tSpp\n');


%gpldatout = gpldat(gpldat.Channel==5,:);
gpldatout = gpldat(k2,:);
writetable(gpldatout(:,[1:9, end]), fname, 'Delimiter', '\t', 'WriteVariableNames', false)    

% export the file
% - Read input file data content (without header).
 fid_in = fopen(fname, 'r') ;             % Open input file for reading.
 fgetl(fid_in) ;                                 % Skip header line in input file.
 content = fread(fid_in) ;                       % Read rest of input file.
 fclose(fid_in) ;                                % Close input file.
 
 
 % - Write output file: new header + previous data content.
 fid_out = fopen(fname, 'w') ;  % Open input file for writing.
 fprintf(fid_out, header') ;       % Output new header.
 fwrite(fid_out, content) ;                      % Output previous data content.
 fclose(fid_out) ;  







