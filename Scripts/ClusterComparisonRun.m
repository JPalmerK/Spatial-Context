
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
localize_struct = trimlocalize_struct(localize_struct, 1, 5)

localize_struct = trimlocalize_struct(localize_struct, 1, 8)


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

%% Error Rates - channel 8 as parent


corr_thresh = linspace(0,1, 10);
Propout =[];
NClustout =[];


% Channel 5
parent =8
close all;


calls_arrivals =struct2table(hyd(parent).detection.calls);



examp = struct();
examp.array_struct = localize_struct.hyd(parent).array;
examp.hydrophone_struct = hydrophone_struct;
examp.randomMiss =0;
examp.s = 12;
examp.child_idx = [1:8];
examp.localize_struct =localize_struct;
examp.limitTime =24*60*60;
examp.maxEltTime = 90;
examp.calls =calls_arrivals;
examp.callParm =  hyd(parent).detection.parm;
examp.fs =2000;
examp.PosUncertsigma = 0.0004^2 +.1^2 + .3^2;
examp.drift = 0;
examp.c =1500;


examp.arrivalTable = updateArrTableGPL(examp); % Run this first!
examp.arrivalArray = UpdateArrArray(examp);
examp.TDOA_vals = UpdateTDOA(examp);
examp.Sim_mat= simMatMaxofProd(examp);



%UpdateTDOA(examp)
%simMatIdealNewSim(examp); % Simulation matrix spatial method

for ii=1:length(corr_thresh)
    sim_thresh=corr_thresh(ii);
    examp.cutoff = sim_thresh;
    examp.chains =updateChainsEncounterFirst(examp);
    examp.Cluster_id= updateClusterID(examp);
    [propCorrect, NClusters] = RunComp(examp, parent, hyd, truth,sim_thresh);
    Propout(ii)=propCorrect;
    NClustout(ii)=NClusters;
    
end



% Do the same for the TDOA only method
%examp.clearCalcValues
%UpdateArrTable(examp) % Run this first!
examp.Sim_mat= simMatTDOAonly(examp);
%examp.simMatTDOAonly()

PropoutTDOA =[];
NClustoutTDOA =[];
for ii=1:length(corr_thresh)
    sim_thresh=corr_thresh(ii);
    examp.cutoff = sim_thresh;
    examp.chains =updateChainsEncounterFirst(examp);
    examp.Cluster_id= updateClusterID(examp);
    [propCorrect, NClusters] = RunComp(examp, parent, hyd, truth,sim_thresh);
    PropoutTDOA(ii)=propCorrect;
    NClustoutTDOA(ii)=NClusters;

end


% Look at TOA only 

examp.Cluster_id = acEnc(examp);
[propCorrectBase, NClusters] = RunComp(examp, parent, hyd, truth,sim_thresh);

figure(3)
subplot(2,1,1)
hold on
plot(NClustoutTDOA, PropoutTDOA, '-+')
plot(NClustoutTDOA, Propout, '-*')
plot(NClustoutTDOA, repmat(propCorrectBase, size(corr_thresh)), '-.k')
xlabel('Number of Acoustic Encounters')
ylabel('Proportion of Calls Correctly Classified')
title('Hydrophone 8')
legend('TDOA Method', 'Spatial Method', 'Baseline')
ylim([.75 1])





%% Error Rates - channel 5 as parent
% Populate data and parameters
parent =5

calls_arrivals =struct2table(hyd(parent).detection.calls);
examp = struct();
examp.array_struct = localize_struct.hyd(parent).array;
examp.hydrophone_struct = hydrophone_struct;
examp.time_cut = 20*60;
examp.randomMiss =0;
examp.s = 12;
examp.child_idx = [1:8];
examp.localize_struct =localize_struct;
examp.limitTime =4*60*60;
examp.maxEltTime = 90;
examp.calls =calls_arrivals;
examp.callParm =  hyd(parent).detection.parm;
examp.fs =2000;
examp.PosUncertsigma = 0.0004^2 +.1^2 + .3^2;
examp.drift = 0;
examp.c =1500;


examp.arrivalTable = updateArrTableGPL(examp); % Run this first!
examp.arrivalArray = UpdateArrArray(examp);
examp.TDOA_vals = UpdateTDOA(examp);

% Spatial Method
examp.Sim_mat= simMatMaxofProd(examp);

for ii=1:length(corr_thresh)
    sim_thresh=corr_thresh(ii);
    examp.cutoff = sim_thresh;
    examp.chains =updateChainsEncounterFirst(examp);
    examp.Cluster_id= updateClusterID(examp);
    [propCorrect, NClusters] = RunComp(examp, parent, hyd, truth,sim_thresh);
    Propout(ii)=propCorrect;
    NClustout(ii)=NClusters;
    
end



% Do the same for the TDOA only method
examp.Sim_mat= simMatTDOAonly(examp);


PropoutTDOA =[];
NClustoutTDOA =[];
for ii=1:length(corr_thresh)
    sim_thresh=corr_thresh(ii);
    examp.cutoff = sim_thresh;
    examp.chains =updateChainsEncounterFirst(examp);
    examp.Cluster_id= updateClusterID(examp);
    [propCorrect, NClusters] = RunComp(examp, parent, hyd, truth,sim_thresh);
    PropoutTDOA(ii)=propCorrect;
    NClustoutTDOA(ii)=NClusters;

end




% TOA/acoustic encounters only
examp.Cluster_id = acEnc(examp);
[propCorrectBase, NClusters] = RunComp(examp, parent, hyd, truth,sim_thresh);



figure(3)
subplot(2,1,2)
hold on
plot(NClustoutTDOA, PropoutTDOA, '-+')
plot(NClustoutTDOA, Propout, '-*')
plot(NClustoutTDOA, repmat(propCorrectBase, size(corr_thresh)), '-.k')
xlabel('Number of Acoustic Encounters')
ylabel('Proportion of Calls Correctly Classified')
title('Hydrophone 5')
legend('TDOA Method', 'Spatial Method', 'Baseline')
ylim([.75 1])





figure(3)
subplot(2,1,2)
hold on
legend('TDOA Method', 'Spatial Method', 'Baseline')
ylim([.75 1])














