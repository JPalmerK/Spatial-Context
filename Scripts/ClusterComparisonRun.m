
% Run experiments
close all; clear all
clear classes; clc

whereAmI = loadSimspace();
dclde_2013_meta = xlsread(whereAmI{1});
load(whereAmI{2})
load(whereAmI{3})
load(whereAmI{4})
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

% Merge the array structure with the localize structure
 
 for i = 1:length(localize_struct.hyd)
    localize_struct.hyd(i).array = array_struct_data(i).array;
 end


% Trim the LSQ space
localize_struct = trimlocalize_struct(localize_struct, .01, 5)

localize_struct = trimlocalize_struct(localize_struct, .01, 8)


%% Pick a parent and go

fs = 2000;
ssp = localize_struct.parm.ssp;
grid_depth = localize_struct.parm.grid_depth;






%% Load truth dataset
% Find the indicies that coinced with the truth table
% Load the calls
%truth = readtable('C:\Users\Kaitlin\Desktop\NOPP6_EST_20090328.selections.txt');
truth = readtable('/home/kpalmer/Desktop/NOPP6_EST_20090328.selections.txt');
truth.mstart = datenum('20090328', 'yyyymmdd')+truth.BeginTime_s_/86400;
truth.mend = datenum('20090328', 'yyyymmdd')+truth.EndTime_s_/86400;
truth.mid = (truth.mend+truth.mstart)/2;


% Diagnostic plots
Chan5data = truth(truth.Channel==5,:);
Chan8data = truth(truth.Channel==8,:);

figure

subplot(2,1,2)
gscatter(Chan8data.BeginTime_s_, zeros(height(Chan8data)),  Chan8data.Species)

subplot(2,1,1)
gscatter(Chan5data.BeginTime_s_, zeros(height(Chan5data)),  Chan5data.Species)

%% Correlate Clusters


% Channel 5
parent =8
close all;


calls_arrivals =struct2table(hyd(parent).detection.calls);


% Populate data and parameters
examp = clusterClassGPL();
examp.array_struct = localize_struct.hyd(parent).array;
examp.hydrophone_struct = hydrophone_struct;
examp.time_cut = 20*60;
examp.randomMiss =0;
examp.s = 12;
examp.child_idx = [1:8];
examp.localize_struct =localize_struct;
examp.limitTime =4*60*60;
examp.maxEltTime =60;
examp.calls =calls_arrivals;
examp.callParm =  hyd(parent).detection.parm;


corr_thresh = linspace(0,1, 20);
Propout =[];
NClustout =[];

calls_arrivals =struct2table(hyd(parent).detection.calls);

examp.clearCalcValues
UpdateArrTable(examp) % Run this first!
examp.simMatIdealNewSim(); % Simulation matrix spatial method

for ii=1:length(corr_thresh)
sim_thresh=corr_thresh(ii);
    [propCorrect, NClusters] = RunComp(examp, parent, hyd, truth,sim_thresh);
    Propout(ii)=propCorrect;
    NClustout(ii)=NClusters;

end



% Do the same for the TDOA only method
examp.clearCalcValues
UpdateArrTable(examp) % Run this first!
examp.simMatTDOAonly()

PropoutTDOA =[];
NClustoutTDOA =[];
for ii=1:length(corr_thresh)
    sim_thresh=corr_thresh(ii);
    [propCorrect, NClusters] = RunComp(examp, parent, hyd, truth,sim_thresh);
    PropoutTDOA(ii)=propCorrect;
    NClustoutTDOA(ii)=NClusters;

end


% Look at TOA only 

examp.clearCalcValues
examp.toaOnlyCluster()
[propCorrectBase, NClusters] = RunComp(examp, parent, hyd, truth,sim_thresh);

figure(2)
subplot(2,1,1)
hold on
plot(corr_thresh, PropoutTDOA, '-+')
plot(corr_thresh, Propout, '-*')
plot(corr_thresh, repmat(propCorrectBase, size(corr_thresh)), '-.k')
xlabel('Similarity Threshold')
ylabel('Proportion of Calls Correctly Classified')
title('Hydrophone 8')
legend('TDOA Method', 'Spatial Method', 'Baseline')






% Populate data and parameters
parent =5

calls_arrivals =struct2table(hyd(parent).detection.calls);
examp = clusterClassGPL();
examp.array_struct = localize_struct.hyd(parent).array;
examp.hydrophone_struct = hydrophone_struct;
examp.time_cut = 20*60;
examp.randomMiss =0;
examp.s = 12;
examp.child_idx = [1:8];
examp.localize_struct =localize_struct;
examp.limitTime =4*60*60;
examp.maxEltTime =60;
examp.calls =calls_arrivals;
examp.callParm =  hyd(parent).detection.parm;
examp.clearCalcValues
UpdateArrTable(examp) % Run this first!
examp.simMatIdealNewSim();

calls_arrivals =struct2table(hyd(parent).detection.calls);

for ii=1:length(corr_thresh)
sim_thresh=corr_thresh(ii);
    [propCorrect, NClusters] = RunComp(examp, parent, hyd, truth,sim_thresh);
    Propout(ii)=propCorrect;
    NClustout(ii)=NClusters;

end



% Do the same for the TDOA only method
examp.clearCalcValues
UpdateArrTable(examp) % Run this first!
examp.simMatTDOAonly()

PropoutTDOA =[];
NClustoutTDOA =[];
for ii=1:length(corr_thresh)
    sim_thresh=corr_thresh(ii);
    [propCorrect, NClusters] = RunComp(examp, parent, hyd, truth,sim_thresh);
    PropoutTDOA(ii)=propCorrect;
    NClustoutTDOA(ii)=NClusters;

end




% Look at TOA only 
examp.clearCalcValues
examp.toaOnlyCluster()

[propCorrectBase, NClusters] = RunComp(examp, parent, hyd, truth,sim_thresh);


figure(2)
subplot(2,1,2)
hold on
plot(corr_thresh, PropoutTDOA, '-+')
plot(corr_thresh, Propout, '-*')
plot(corr_thresh, repmat(propCorrectBase, size(corr_thresh)), '-.k')
xlabel('Similarity Threshold')
ylabel('Proportion of Calls Correctly Classified')
title('Hydrophone 5')
legend('TDOA Method', 'Spatial Method', 'Baseline')





figure(2)
subplot(2,1,2)
hold on
legend('TDOA Method', 'Spatial Method', 'Baseline')














