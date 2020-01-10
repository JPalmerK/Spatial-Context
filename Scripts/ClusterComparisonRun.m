
% Run experiments
close all; clear all; clc

% Load metadata (deployment information), GPL localize structure (has TDOA
% values etc.) hyd structure (origional detections, dex in localize
% strucutre references call rows in hyd) and array_struct_data which has
% array structures created with hydrophoens 5 and 8 run as the 'parent' or
% 'master' by GPL
dclde_2013_meta = xlsread('DCL2013_NEFSC_SBNMS_200903_metadata.xlsx');
load('DCLDE2013_DCLDE_2013_10_Chan201914101_96_localize_struct.mat')
load('DCLDE2013_RW_DCLDE_2013_10_Chan201908101_96_hyd.mat')
load('DCLDE2013_RWDCLDE_2013_10_Chan201914101_96_array_struct_data.mat')
clear whereAmI

% Convert meta data to what GPL expects
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

localize_struct.hyd(5).detectorScore= zeros(...
    [1 length(localize_struct.hyd(5).dex)])/0;

parents = [5,8];

% Relate the detector score to the localize structure
for jj =1:length(parents)
    parent = parents(jj);
    for ii=1:length(localize_struct.hyd(parent).dex)
        
        score = hyd(parent).detection.calls(localize_struct.hyd(parent).dex(ii)).calltype_1_score;1
        localize_struct.hyd(parent).detectorScore(ii)=score;
    end
end

% Pre-compute ambiguity surfaces was run with the following parameters and
% as such, more lienient parameters will not work
% Template correlation threshold: 0.01
% Correaltion threhsold (correlation between channels): 0.25
% MAX TDOA (trimming TDOA where values are outside of the expected range)
% maximum expected TDOA plus 0.5 sec clock drift
% 18km maximum detection range (must be detected by two sensors, swim speed
% 6m/s


% Load truth dataset
% Find the indicies that coinced with the truth table
% Load the calls
truth = readtable('NOPP6_EST_20090328.selections.txt');
%truth = readtable('/home/kpalmer/Desktop/NOPP6_EST_20090328.selections.txt');
truth.mstart = datenum('20090328', 'yyyymmdd')+truth.BeginTime_s_/86400;
truth.mend = datenum('20090328', 'yyyymmdd')+truth.EndTime_s_/86400;
truth.mid = (truth.mend+truth.mstart)/2;
truth.Sec = (truth.mstart-min(truth.mstart))*60*60*24;

% Link validation spreadsheet to the detector output
for ii=1:length(parents)
    parent = parents(ii);
    
    % Create a validation column (T.H. code calls this 'pruned' sticking
    % with the naming scheme)
    localize_struct.hyd(parent).pruned = ...
        zeros(size(localize_struct.hyd(parent).dex))/0;
    
    % Create temporary subset of the truth labels
    truthSub = truth(truth.Channel==parent,:);
    
    detected=0;
    % Step through the calls with TDOA info and if they match a validation
    % call assign 1 to the pruned value if it's a rw and 0 otherwise.
    for jj=1:length(localize_struct.hyd(parents(ii)).dex)
        
        % To get the call start get the index(dex) from the localize struct
        % and reference it back to thy hyd file (detections)
        callDex = localize_struct.hyd(parent).dex(jj);
        
        % Call start time in Matlab form (<3 Matlab)
        callStart = hyd(parent).detection.calls(callDex).julian_start_time;
        
        % See if any calls are within .3 of a seconds from a truth label
        [diffval, diffIdx] = min(abs(truthSub.mstart - callStart));
        diffSec = diffval*24*60*60;
        

        
        % If a detection correlates with the truth value add the pruned fx
        if diffSec<.5
            %disp('blarg')
            
            % Species id
            spp = truthSub.Species(diffIdx);
            
            if strcmp(spp, 'rw')
                localize_struct.hyd(parent).pruned(jj) = 1;
            else
                localize_struct.hyd(parent).pruned(jj) = 0;
            end
            truthSub(diffIdx,:)=[]; 
            
        end
    end

    clear truthSub
end


% Determine how many raw detections coincided with rw annotations
for ii=1:length(parents)
    detected=0;
    parent = parents(ii);
    
    % Create temporary subset of the truth labels
    truthSub = truth(truth.Channel==parent,:);
    callStart = cell2mat({hyd(parent).detection.calls.julian_start_time});
    for jj=1:height(truthSub)
        % Differences from the origional detection (for making results table)
        [diffvaltmp, diffIdxtmp] = min(abs(truthSub.mstart(jj)-...
            callStart));
        
        if (diffvaltmp*24*60*60)<.5
            spp = truthSub.Species(jj);
            if strcmp(spp, 'rw')
                detected = detected+1;
            end
        end
    end
    
end

% Create subset of truth data where calls are interleaved
% truth.Sec = (truth.mstart-floor(truth.mstart))*60*60*24;
% truth_sub = truth(truth.BeginTime_s_>=6.49e04 &...
%     truth.BeginTime_s_<= 6.85e04,:);
% 
% startCut = truth.BeginTime_s_(truth.Selection==3045);
% arr_times_sec = localize_struct.hyd(5).rtimes/2000;
% noninterleavedIdx = find(arr_times_sec<startCut);
% 
% localize_struct.hyd(5).delays(noninterleavedIdx,:)=[];
% localize_struct.hyd(5).cross_score(noninterleavedIdx,:)=[];
% localize_struct.hyd(5).coord_time(noninterleavedIdx,:)=[];
% localize_struct.hyd(5).rtimes(noninterleavedIdx)=[];
% localize_struct.hyd(5).dex(noninterleavedIdx)=[];
% localize_struct.hyd(5).coordinates(:,:,noninterleavedIdx)=[];
% localize_struct.hyd(5).score(:,noninterleavedIdx)=[];
% localize_struct.hyd(5).detectorScore(noninterleavedIdx)=[];
% 
% localize_struct.hyd(5).pruned(noninterleavedIdx)=[];



% Histograms of the scores for the TP and FP detections
figure
labelSpp=[{'Noise'}, {'Right Whale'}]
idxvals =[0 1]
for ii=1:length(idxvals)
    
    idx = find(localize_struct.hyd(5).pruned == idxvals(ii));
    scores = localize_struct.hyd(5).detectorScore(idx);
    
    subplot(1,2,ii); hist((scores),[0:.05:1]);
    title([labelSpp(ii)]);
    xlabel('Template Detection Scores');
    median(scores)
end



%% Do some trimming of the localization structure to remove erronious detections and associations

% Trim the localization structure based on the template correlation
localize_struct_trimmed = trimlocalize_struct(localize_struct,...
    hyd, .01);


% Trim the localizeation structure based on the call association score
 localize_struct_trimmed = trimlocalize_structCrossCorr(...
     localize_struct_trimmed, hyd, .25);


% % % Remove all the delays where the TDOA is greater than the array geometry
localize_struct_trimmed = trimlocalize_structTDOA(localize_struct_trimmed,...
    hyd, 0);



% % Remove rows from the localizatoin structure where the calls are not
% % detected on two or more hydrophones
% Hydrophone 5 as parent - idxs= 1,2,3,4,6

badidx = find(all...
    (isnan(localize_struct_trimmed.hyd(5).delays(:, [1,2,3,4,6])),2));


localize_struct_trimmed.hyd(5).delays(badidx,:)=[];
localize_struct_trimmed.hyd(5).cross_score(badidx,:)=[];
localize_struct_trimmed.hyd(5).coord_time(badidx,:)=[];
localize_struct_trimmed.hyd(5).rtimes(badidx)=[];
localize_struct_trimmed.hyd(5).dex(badidx)=[];
localize_struct_trimmed.hyd(5).coordinates(:,:,badidx)=[];
localize_struct_trimmed.hyd(5).score(:,badidx)=[];
localize_struct_trimmed.hyd(5).detectorScore(badidx)=[];
%localize_struct_trimmed.hyd(5).pruned(badidx)=[];
%%



% Remove indicies where no call associated
nanIdx = find(isnan(localize_struct_trimmed.hyd(5).pruned));
localize_struct_trimmed.hyd(5).delays(nanIdx,:)=[];
localize_struct_trimmed.hyd(5).cross_score(nanIdx,:)=[];
localize_struct_trimmed.hyd(5).coord_time(nanIdx,:)=[];
localize_struct_trimmed.hyd(5).rtimes(nanIdx)=[];
localize_struct_trimmed.hyd(5).dex(nanIdx)=[];
localize_struct_trimmed.hyd(5).coordinates(:,:,nanIdx)=[];
localize_struct_trimmed.hyd(5).score(:,nanIdx)=[];
localize_struct_trimmed.hyd(5).detectorScore(nanIdx)=[];

localize_struct_trimmed.hyd(5).pruned(nanIdx)=[];


% Get the total number of detections and TP for precision/recall

figure
labelSpp=[{'False Positive'}, {'Right Whale'}]
idxvals =[0 1 ]
for ii=1:length(idxvals)
    
    idx = find(localize_struct_trimmed.hyd.pruned == idxvals(ii));
    scores = localize_struct_trimmed.hyd.detectorScore(idx);
    
    subplot(1,3,ii); hist((scores),[0:.05:1]);
    title([labelSpp(ii)]);
    xlabel('Template Detection Scores');
end


%% Create structure (formally class) for running the analysis

fs = 2000;
ssp = localize_struct.parm.ssp;
grid_depth = localize_struct.parm.grid_depth;

parent = 5;
examp = struct();
examp.hydrophone_struct = hydrophone_struct;
examp.randomMiss =0;
examp.s = 6; % maximum swim speed of the animal in m/s
examp.child_idx = [1,2,3,4,6]; % local array for parent 5, [1,2,3,4,6], for 8 master 4:8
examp.fs =fs;
examp.PosUncertsigma = 0.0004^2 +.1^2 + .3^2; % seconds see EM Nosal 20
examp.drift = .5; % clock drift uncertainty
examp.c =1500;  % speed of sound
examp.truncateKm=18; % replaced by half normal function, previously maximum detection distance
examp.array_struct = array_struct_data(parent).array;

% where the projected ambiguity surfaces are stored or will be stored if
% the file location is emptyish
examp.propAmbLoc = 'D:\DCLDE2013_ProjAmbSurfs'; 

% Create filter grids to knock out ambiguity surfaces where calls can't
% have come from (half normal function)
examp.filtGrid = createDetRangeFiltGrid(examp, hydrophone_struct);

% try: 
% figure; imagesc(examp.array_struct.latgrid, ...
% examp.array_struct.latgrid, examp.filtGrid(:,:,1)), axis xy




%% Create the structures to house the information (formally objects, boo!)

examp.localize_struct =localize_struct_trimmed;
examp.maxEltTime = 20*60; % Maximum elapsed time between detections to make another encounter
examp.callParm =  hyd(parent).detection.parm; % From GPL
examp.arrivalTable = updateArrTableGPL(examp); % Run this first!
examp.arrivalArray = UpdateArrArrayRealData(examp);
examp.TDOA_vals = UpdateTDOA(examp);

% Create coppies for the different method comparisons
exampSpatial = examp;
exampTDOA = examp;
exampBaseline = examp;
Results=struct();

%% Create simulation matricies (extra large, then trim down as time threshs)

% Create the similarity matrix using pre-computed ambiguity surfaces
simMatTemp = simMatMaxofProdPreComputed(examp);
exampSpatial.Sim_mat= simMatTemp;

% Create the similarity matrix using the TDOA approach
exampTDOA.Sim_mat = simMatTDOAonly(examp);

% Create the similarity matrix (it's all 1's) using the acoustic encounters
% only approach
exampBaseline.Sim_mat = ones(size(exampTDOA.Sim_mat));
exampBaseline.Cluster_id = acEnc(examp);
%% Create the results field to look at precision and recall as a function 
% of the max elapsed time and similarity threshold

simThreshs = linspace(.1,.98, 9);
TimeThreshs = fliplr(linspace(8,120,10));

ExperimentTDOA = struct();
ExperimentMaxProd = struct();
ExperimentUnaided = struct();

% Pre allocate the output matrix
PrecisionMatSpatial = inf(length(simThreshs), length(TimeThreshs),...
    length(simThreshs));
RecallMatSpatial = PrecisionMatSpatial;
PrecisionMatTDOA = PrecisionMatSpatial;
RecallMatTDOA = PrecisionMatSpatial;
PrecisionMatBaseline =PrecisionMatSpatial;
RecallMatBaseline =PrecisionMatSpatial;
PrecisionMatAcEnc= PrecisionMatSpatial;
RecallMatAcEnc =PrecisionMatSpatial;
% Create the basic structure, this needs to be done for every correlation
% threshold

nmi_baseline = zeros(length(simThreshs), length(TimeThreshs));
nmi_tdoa = nmi_baseline;
nmi_spatial = nmi_baseline;

for ii = 1:length(simThreshs)
    simThresh = simThreshs(ii);
    exampSpatial.cutoff = simThresh;
    exampTDOA.cutoff = simThresh;
    exampBaseline.cutoff=.5;
    disp('next simThresh')
    
    for jj = 1:length(TimeThreshs)
        % Create precision recall curves for each similarity and time threshold
        
            exampSpatial.maxEltTime = TimeThreshs(jj);
            exampTDOA.maxEltTime = TimeThreshs(jj);
            exampBaseline.maxEltTime = TimeThreshs(jj);
           
            
            % Run the clustering spatial
            exampSpatial.chains = updateChainsEncounterFirst(exampSpatial);
            exampSpatial.Cluster_id= updateClusterID(exampSpatial);
            
            %Run the clustering TDOA
            exampTDOA.chains = updateChainsEncounterFirst(exampTDOA);
            exampTDOA.Cluster_id= updateClusterID(exampTDOA);
            
            % Baseline clustering, time only
            exampBaseline.chains = updateChainsEncounterFirst(exampBaseline);
            exampBaseline.Cluster_id= updateClusterID(exampBaseline);
            
            % Spatial and baseline
            [ClusteredSpatial, Detector,...
                gpldatout_chan, nmi ] = ...
                PrecisionRecallGPL(exampSpatial, parent, hyd);
            Results.Spatial(ii,jj) = ClusteredSpatial;
            Results.SpatialClusters(ii,jj) = {exampSpatial.Cluster_id};
            Results.GPLTable(ii,jj) = {gpldatout_chan};
            nmi_spatial(ii,jj) =nmi;
            
            if ii == jj ==1
                Results.Baseline(1,1) = Detector;
            end
            
            % TDOA precision/recall
            [ClusteredTDOA, ~,~,nmi] = PrecisionRecallGPL(exampTDOA, parent,...
                hyd);
            Results.TDOA(ii,jj) = ClusteredTDOA;
            Results.TDOAClusters(ii,jj) = {exampTDOA.Cluster_id};
            nmi_tdoa(ii,jj) =nmi;
            
            
            % Baseline clustering precision recall
            [ClusteredBase, ~,~,nmi] = PrecisionRecallGPL(exampBaseline, parent,...
                hyd);
            Results.AcousticEncounters(ii,jj) = ClusteredBase;
            Results.AcousticEncountersClusters(ii,jj) = ...
                {exampBaseline.Cluster_id};
            nmi_baseline(ii,jj) =nmi;
            disp('next correlatinthresh')
        
    end

    disp('Finish Sim Thres Iter')
end

%% Create a nice plot

simIdx =2;
simIdx2=5;
simIdx3=9;

figure
hold on
Markers = {'-o','-.','-+','-*'};
sizeVal = [2, 6,4,2];
markercounter=1;


for jj=1:3:length(TimeThreshs)
    timeidx =jj;
    

        subplot(3,3,1)
        hold on
        plot(Results.TDOA(simIdx,timeidx).Recall,...
            Results.TDOA(simIdx,timeidx).Precision, ...
            Markers{markercounter},'color',[230/255 159/255 0/255],...
            'MarkerSize',sizeVal(markercounter))
        if jj==10
         plot(Results.Baseline(1,1).Recall,...
            Results.Baseline(1,1).Precision, 'k-','Linewidth',2)
        end
        ylabel('Precision');
        title(['TDOA: ', num2str(round(simThreshs(simIdx),2)), ' Sim. Thresh'])
        
        subplot(3,3,2)
        hold on
        plot(Results.Spatial(simIdx,timeidx).Recall,...
            Results.Spatial(simIdx,timeidx).Precision,...
            Markers{markercounter},'color',[86/255 180/255 233/255],...
            'MarkerSize',sizeVal(markercounter))
        if jj==10
         plot(Results.Baseline(1,1).Recall,...
            Results.Baseline(1,1).Precision, 'k-','Linewidth',2)
        end
        title(['Ambiguity Surfaces: ',...
            num2str(round(simThreshs(simIdx),2)), ' Sim. Thresh'])
        
        subplot(3,3,3)
        hold on
        plot(Results.AcousticEncounters(simIdx,timeidx).Recall,...
            Results.AcousticEncounters(simIdx,timeidx).Precision,...
             Markers{markercounter},'color',[0/255 158/255 115/255],...
            'MarkerSize',sizeVal(markercounter))
        if jj==10
         plot(Results.Baseline(1,1).Recall,...
            Results.Baseline(1,1).Precision, 'k-','Linewidth',2)
        end
        title('Acoustic Encounters')
        
        
        % Second Row
        
        subplot(3,3,4)
        hold on
        plot(Results.TDOA(simIdx2,timeidx).Recall,...
            Results.TDOA(simIdx2,timeidx).Precision, ...
            Markers{markercounter},'color',[230/255 159/255 0/255],...
            'MarkerSize',sizeVal(markercounter))
        if jj==10
         plot(Results.Baseline(1,1).Recall,...
            Results.Baseline(1,1).Precision, 'k-','Linewidth',2)
        end
        ylabel('Precision');
        title(['TDOA: ',num2str(round(simThreshs(simIdx2),2)), ' Sim. Thresh'])
        
        subplot(3,3,5)
        hold on
        plot(Results.Spatial(simIdx2,timeidx).Recall,...
            Results.Spatial(simIdx2,timeidx).Precision,...
            Markers{markercounter},'color',[86/255 180/255 233/255],...
            'MarkerSize',sizeVal(markercounter))
        if jj==10
         plot(Results.Baseline(1,1).Recall,...
            Results.Baseline(1,1).Precision, 'k-','Linewidth',2)
        end
        title(['Ambiguity Surfaces: ',num2str(round(simThreshs(simIdx2),2)), ' Sim. Thresh'])
        
        subplot(3,3,6)
        hold on
        plot(Results.AcousticEncounters(simIdx2,timeidx).Recall,...
            Results.AcousticEncounters(simIdx2,timeidx).Precision,...
             Markers{markercounter},'color',[0/255 158/255 115/255],...
            'MarkerSize',sizeVal(markercounter))
        if jj==10
         plot(Results.Baseline(1,1).Recall,...
            Results.Baseline(1,1).Precision, 'k-','Linewidth',2)
        end
        title('Acoustic Encounters')
        
         
        % Third Row
        subplot(3,3,7)
           hold on
        plot(Results.TDOA(simIdx3,timeidx).Recall,...
            Results.TDOA(simIdx3,timeidx).Precision, ...
            Markers{markercounter},'color',[230/255 159/255 0/255],...
            'MarkerSize',sizeVal(markercounter))
        if jj==10
         plot(Results.Baseline(1,1).Recall,...
            Results.Baseline(1,1).Precision, 'k-','Linewidth',2)
        end
         xlabel('Recall'); ylabel('Precision');
        title(['TDOA: ', num2str(round(simThreshs(simIdx3),2)), ' Sim. Thresh'])
        
        subplot(3,3,8)
        hold on
        plot(Results.Spatial(simIdx3,timeidx).Recall,...
            Results.Spatial(simIdx3,timeidx).Precision,...
            Markers{markercounter},'color',[86/255 180/255 233/255],...
            'MarkerSize',sizeVal(markercounter))
        if jj==10
         plot(Results.Baseline(1,1).Recall,...
            Results.Baseline(1,1).Precision, 'k-','Linewidth',2)
        end
         xlabel('Recall');
        title(['Ambiguity Surfaces: ', num2str(round(simThreshs(simIdx3),2)), ' Sim. Thresh'])
        
        subplot(3,3,9)
        hold on
        plot(Results.AcousticEncounters(simIdx3,timeidx).Recall,...
            Results.AcousticEncounters(simIdx3,timeidx).Precision,...
             Markers{markercounter},'color',[0/255 158/255 115/255],...
            'MarkerSize',sizeVal(markercounter))
        if jj==10
         plot(Results.Baseline(1,1).Recall,...
            Results.Baseline(1,1).Precision, 'k-','Linewidth',2)
        end
        xlabel('Recall');
        title('Acoustic Encounters')
        
              
        markercounter = markercounter+1;
        
         
         
    
end


legendtxt =[{'120 s'}, {'82 s'}, {'45 s'}, {'8 s'}, {'Template Corr.'}];
% Add the legends
for ii=1:9
    subplot(3,3,ii)
    legend(legendtxt, 'Location','southwest')
end


%%
close all
jitterAmount = 0.000;

for jj=1:length(TimeThreshs)
    timeidx =jj;
    figure;
    for ii =1:length(simThreshs)
        
        titlestr =[num2str(TimeThreshs(timeidx),'%0.2f'), ' sec TimeThresh ',...
            num2str(num2str(simThreshs(ii),'%0.2f')), ' SimThresh'];
        
        subplot(3,3,ii)
        
        hold on
        
        % Jitter the baseline detector so it shows up
        jitterValuesX = 2*(rand(size(Results.Baseline(1,1)))-0.5)...
            *jitterAmount;   % +/-jitterAmount max
        jitterValuesY = 2*(rand(size(Results.Baseline(1,1).Precision))-0.5)...
            *jitterAmount;   % +/-jitterAmount max
        
        
        plot(Results.TDOA(ii,timeidx).Recall,...
            Results.TDOA(ii,timeidx).Precision, '.-')
        plot(Results.Spatial(ii,timeidx).Recall,...
            Results.Spatial(ii,timeidx).Precision, '.-')
        plot(Results.AcousticEncounters(ii,timeidx).Recall,...
            Results.AcousticEncounters(ii,timeidx).Precision, '.-')
        
        plot(Results.Baseline(1,1).Recall+jitterValuesX,...
            Results.Baseline(1,1).Precision+jitterValuesY, '.-')
        xlabel('Recall'); ylabel('Precision');
        
        title(titlestr)
        
        
    end
end

%% Make plots showing cluster accuracy

 %legend('TDOA', 'Spatial', 'Acoustic Encounters Only','Detector Only')
 
 
figure
subplot(2,2,1); imagesc(TimeThreshs, simThreshs,nmi_baseline); colorbar
title('Acoustic Encounters')
xlabel('Elapsed Time Threshold (s)'); ylabel('Similarity Threshold');
%caxis([.2 .27])

subplot(2,2,2); imagesc(TimeThreshs, simThreshs,nmi_tdoa); colorbar
title('TDOA')
xlabel('Elapsed Time Threshold (s)'); ylabel('Similarity Threshold');
%caxis([.2 .27])

subplot(2,2,3); imagesc(TimeThreshs, simThreshs,nmi_spatial); colorbar
title('Spatial')
xlabel('Elapsed Time Threshold (s)'); ylabel('Similarity Threshold');
%caxis([.2 .27]) 
 
 

%% Area under the precision recall curve
maxA1 =0; maxA2 =0; maxA3=0;
figure
simIdx = 1;
for ii=1:length(TimeThreshs)

    subplot(1,3,1)
    hold on
    area =Results.TDOA(simIdx,ii).Recall.*...
        Results.TDOA(simIdx,ii ).Precision;
     maxA1 = max([area, maxA1]);
     plot(area, '.-')
     plot(Results.Baseline.Recall .*Results.Baseline.Precision, 'k')
    title('TDOA')
    ylabel('Area under the Precision/Recall Curve')
    
    subplot(1,3,2)
    hold on
    
    area =Results.Spatial(simIdx,ii).Recall.*...
        Results.Spatial(simIdx,ii ).Precision;
    maxA2 = max([area, maxA2]);
    plot(area, '.-')
         plot(Results.Baseline.Recall .*Results.Baseline.Precision, 'k')
    %ylim([.4 .8])
    
    title('Spatial')
    
    subplot(1,3,3)
    hold on
    
    area =Results.AcousticEncounters(simIdx,ii).Recall.*...
        Results.AcousticEncounters(simIdx,ii ).Precision;
     maxA3 = max([area, maxA3]);
          plot(Results.Baseline.Recall .*Results.Baseline.Precision, 'k')
     plot(area, '.-')
    %ylim([.4 .8])
    
    title('Acoustic Encounters')
    
    maxOut =[maxA1 maxA2 maxA3];

end


