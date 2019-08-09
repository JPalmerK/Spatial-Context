% Run experiments
close all; clear all; clc

% Load the workspace (not KJP)
load('NewCompMethod.mat')

% % Load the workspace(KJP)
% whereAmI = loadSimspace();
% dclde_2013_meta = xlsread(whereAmI{1});
% load(whereAmI{2});

grid1= cell2mat( array_struct.toa_diff(2));
grid2 =  cell2mat( array_struct.toa_diff(3));
% **** EVA **** I don't understand these structures, but will proceed
% assuming they are correct

% Get distance in meters between the lower and upper right (deg lat)
grid_v = vdist(min(array_struct.latgrid),...
    min(array_struct.longrid),...
    max(array_struct.latgrid),...
    min(array_struct.longrid));

% Get distance in meters between the lower left and lower right (deg lon)
grid_h = vdist(min(array_struct.latgrid),...
    min(array_struct.longrid),...
    min(array_struct.latgrid),...
    max(array_struct.longrid));

% **** EVA **** What are the units here? m? (With the current data this
% gives a grid spacing of 111 m for lat and ~83m for lon).
deltalat_space = grid_v/ (length(array_struct.latgrid)-1);
deltalon_space = grid_h/ (length(array_struct.longrid)-1);

% *** EVA *** need grid labels in m.
% Doing this temporarily - slighly off due to curvature effects etc but
% good enough for now
grid_xx = (0:length(array_struct.longrid)-1)*deltalon_space;
grid_yy = (0:length(array_struct.latgrid)-1)*deltalat_space;



% How many grid squares per second can the whale move
s =6; % max animal swim speed (m/s)
lat_persec = s / deltalat_space; % Number of grid cells agent can move in one second
lon_persec = s / deltalon_space;

% Maximum change in tdoa per second 
maxDTDOA = sqrt(lat_persec^2+lon_persec^2);

dx = grid_xx(1);
dy = grid_yy(1);

%% Example 1: two calls both detected by two hydrophone pairs

close all 
tdoaObs = [0 0]; % call 1 tdoa observations
gridAdj1 = grid1-tdoaObs(1);
gridAdj2 = grid2-tdoaObs(2);

% Combine on third axis
gridCall1 = cat(3,gridAdj1, gridAdj2);
gridCall1 =  normpdf(gridCall1, 0, sigma).*sigma*sqrt(2*pi); % EVA modified;

% Nan's kluge it up in the real [messy] data
gridCall1(isnan(gridCall1))=.001;
AmbCall1 = prod(gridCall1,3);

% Delta TDOA, simulate second agent moving
times = linspace(0, 60, 15) % simulate 20 points over 10 min

deltatdoa1 = times*(.5* maxDTDOA); % change in TDOA when the agent is moving half speed

score =[]; % similarity scores
dists=[]; % distance travelled


% Iterate through time, move the agent and calculate the similarity score

for ii =1:length(deltatdoa1)
    % second call, shifted a bit
    gridAdj1 = grid1-(tdoaObs(1)+deltatdoa1(ii));
    gridAdj2 = grid2- (tdoaObs(2)+deltatdoa1(ii));
    gridCall1 = cat(3,gridAdj1, gridAdj2);
    gridCall1 =  normpdf(gridCall1, 0, sigma).*sigma*sqrt(2*pi);
    
    % Ambiguity surface call 2
    gridCall1(isnan(gridCall1))=0.01
    AmbCall2 = prod(gridCall1,3);
    
    
    % Similarity score and distance travelled
    [dist, corrScore]= crossCorrSimScores(AmbCall2, AmbCall1, ...
        deltalat_space, deltalon_space);
    
    figure; 
    subplot(1,2,1);
    imagesc(AmbCall1), axis xy, colorbar
    xlabel('relative x-position (km)'), ylabel('relative y-position (km)')
    subplot(1,2,2);
    imagesc(AmbCall2), axis xy, colorbar;
    xlabel('relative x-position (km)'), ylabel('relative y-position (km)')
    score =[score, corrScore];
    dists=[dists, dist];
end
figure;
plot(deltatdoa1, score)
xlabel('delta TDOA')
ylabel('Similarity Score')
title('Two Hyd, Two Hyd- Moving Animal')

%% Example 1: Time Projection only
close all
timegaps = linspace(0,100,10);
tdoaObs = [0 3];
gridAdj1 = grid1-tdoaObs(1);
gridAdj2 = grid2-tdoaObs(2);
gridCall1 = cat(3,gridAdj1, gridAdj2);

% Ambiguity Surface for call 1 (projected in time)
gridCall1 =  normpdf(gridCall1, 0, sigma).*sigma*sqrt(2*pi);

gridCall1(isnan(gridCall1))=.001;
AmbCall1 = prod(gridCall1,3);

% Ambiguity surface for call two that will not be projected in time
AmbCall2 = AmbCall1;

% Size of the filters 
filt_size_lat = ceil(lat_persec * timegaps);
filt_size_lon = ceil(lon_persec * timegaps);
filt_size = ([filt_size_lat; filt_size_lon])';
score=[];

for jj= 1:length(timegaps)
    
    
    % Grow the likelihood space based using image
    % processing max filter. Set the filter size based
    % on the maximum swim speed
    
    AmbCallProj = imdilate(AmbCall1, true(filt_size(jj,:)));
    % Get the comparison value of the projected prob
    % loc space and the next call in the squence
    
    
    
    % Similarity score and distance travelled
    [dist, corrScore]= crossCorrSimScores(AmbCall2, AmbCall1, ...
        deltalat_space, deltalon_space);
    

        figure; 
    subplot(1,2,1);
    imagesc(AmbCall1), axis xy, colorbar
    xlabel('relative x-position (km)'), ylabel('relative y-position (km)')
    subplot(1,2,2);
    imagesc(AmbCallProj), axis xy, colorbar;
    xlabel('relative x-position (km)'), ylabel('relative y-position (km)')
    score =[score, corrScore];
    dists=[dists, dist];
end


figure;
plot(timegaps, score)
xlabel('delta TDOA')
ylabel('Similarity Score')
title('Two Hyd, Two Hyd- Time Proj')



%% Example 2- One call detected 2hyd pairs, one cal detected 1 hyd pair

tdoaObs = [0 3]; % call 1 tdoa observations
gridAdj1 = grid1-tdoaObs(1);
gridAdj2 = grid2-tdoaObs(2);

gridCall1 = cat(3,gridAdj1, gridAdj2);
gridCall1 =  normpdf(gridCall1, 0, sigma).*sigma*sqrt(2*pi);

gridCall1(isnan(gridCall1))=.001;
AmbCall1 = prod(gridCall1,3);


% Delta TDOA, simulate second agent moving
times = linspace(0, 60, 15) % simulate 20 points over 10 min

deltatdoa1 = times*(.5* maxDTDOA); % change in TDOA when the agent is moving half speed


score =[]; % similarity scores
dists=[]; % distance travelled

for ii =1:length(times)
    
    % second call, shifted a bit
    tdoaObs =[tdoaObs(1)+deltatdoa1(ii)];
    
    
    
    % Call 2 grid 1
    AmbCall2 =  normpdf(grid1-tdoaObs, 0, sigma).*sigma*sqrt(2*pi);
    
    
    figure; subplot(1,2,1);imagesc(AmbCall1), axis xy, colorbar
    xlabel('relative x-position (km)'), ylabel('relative y-position (km)')
    subplot(1,2,2);
    imagesc(AmbCall2), axis xy, colorbar
    xlabel('relative x-position (km)'), ylabel('relative y-position (km)')
    
    % Similarity score and distance travelled
    [dist, corrScore]= crossCorrSimScores(AmbCall2, AmbCall1, ...
        deltalat_space, deltalon_space);
    
    score = [score, corrScore];
    dists=[dists, dist];
end

figure; plot(times, score)
xlabel('delta TDOA')
ylabel('Similarity Score')
title('Two Hyd, One Hyd')

%% Example 2- Time projection only
close all
timegaps = linspace(0,100,10);
tdoaObs = [0 3];
gridAdj1 = grid1-tdoaObs(1);
gridAdj2 = grid2-tdoaObs(2);
gridCall1 = cat(3,gridAdj1, gridAdj2);

% Ambiguity Surface for call 1 (projected in time)
gridCall1 =  normpdf(gridCall1, 0, sigma).*sigma*sqrt(2*pi);

gridCall1(isnan(gridCall1))=.001;
AmbCall1 = prod(gridCall1,3);

% Ambiguity surface for call two that will not be projected in time
AmbCall2 = squeeze(gridCall1(:,:,1));



% Size of the filters 
filt_size_lat = ceil(lat_persec * timegaps);
filt_size_lon = ceil(lon_persec * timegaps);
filt_size = ([filt_size_lat; filt_size_lon])';
score=[];

for jj= 1:length(timegaps)
    
   
    % Grow the likelihood space based using image
    % processing max filter. Set the filter size based
    % on the maximum swim speed
    
    AmbCallProj = imdilate(AmbCall1, true(filt_size(jj,:)));
    % Get the comparison value of the projected prob
    % loc space and the next call in the squence
    
  
    figure; subplot(1,2,2);
    imagesc(AmbCallProj), axis xy, colorbar
    xlabel('relative x-position (km)'), ylabel('relative y-position (km)');
    subplot(1,2,1);
    imagesc(AmbCall2), axis xy, colorbar
    xlabel('relative x-position (km)'), ylabel('relative y-position (km)');  
    
    % Similarity score and distance travelled
    [dist, corrScore]= crossCorrSimScores(AmbCall2, AmbCall1, ...
        deltalat_space, deltalon_space);
    
        score =[score, corrScore];
        dists=[dists, dist];
    
%     Second MR method
%     [corrScore]= crossCorrSimScoresFilt(AmbCall2, AmbCallProj, ...
%         filt_size(jj,1), filt_size(jj,2));
%     
%     score =[score, corrScore];
%     dists=[dists, dist];
    
end


figure;
plot(timegaps, score)
xlabel('delta TDOA')
ylabel('Similarity Score')
title('Two Hyd, Two Hyd')


%% Example 3- Two calls, each detected by one hydrophone pair
close all

tdoaObs = [0 3]; % call 1 tdoa observations
gridAdj1 = grid1-tdoaObs(1);
AmbCall1 =  normpdf(gridAdj1, 0, sigma).*sigma*sqrt(2*pi);



% Delta TDOA, simulate second agent moving
times = linspace(0, 60, 15) % simulate 20 points over 10 min

deltatdoa1 = times*(.5* maxDTDOA); % change in TDOA when the agent is moving half speed

score =[]; % similarity scores
dists=[]; % distance travelled

for ii =1:length(deltatdoa1)
    % second call, shifted a bit (moving animal)
    tdoaObsadj =tdoaObs(2)+deltatdoa1(ii);
    AmbCall2 = normpdf(grid2-tdoaObsadj, 0, sigma).*sigma*sqrt(2*pi);

    figure;
    subplot(1,2,1);
    imagesc(AmbCall1), axis xy, colorbar
    xlabel('relative x-position (km)'), ylabel('relative y-position (km)');
    subplot(1,2,2);
    imagesc(AmbCall2), axis xy, colorbar
    xlabel('relative x-position (km)'), ylabel('relative y-position (km)');
    [dist, corrScore]= crossCorrSimScores(AmbCall2,AmbCall1,...
        deltalat_space, deltalon_space);
    score = [score, corrScore];
    dists=[dists, dist];
end

figure
plot(deltatdoa1, score)
xlabel('delta TDOA')
ylabel('Similarity Score')
title('One Hyd, One Hyd')

%% Example 3- Time projection only
close all
timegaps = linspace(0,100,10);
tdoaObs = [0 3];
gridAdj1 = grid1-tdoaObs(1);

% Ambiguity Surface for call 1 (projected in time)
AmbCall1 =  normpdf(gridAdj1, 0, sigma).*sigma*sqrt(2*pi);


% Ambiguity surface for call two that will not be projected in time
AmbCall2 =  normpdf(grid2-tdoaObs(2), 0, sigma).*sigma*sqrt(2*pi);

% Size of the filters 
filt_size_lat = ceil(lat_persec * timegaps);
filt_size_lon = ceil(lon_persec * timegaps);
filt_size = ([filt_size_lat; filt_size_lon])';
score=[];

for jj= 1:length(timegaps)
    
   
    % Grow the likelihood space based using image
    % processing max filter. Set the filter size based
    % on the maximum swim speed 
    
    AmbCallProj = imdilate(AmbCall1, true(filt_size(jj,:)));
    % Get the comparison value of the projected prob
    % loc space and the next call in the squence
    
  
    figure; subplot(1,2,1);
    imagesc(AmbCallProj), axis xy, colorbar
    xlabel('relative x-position (km)'), ylabel('relative y-position (km)');
    subplot(1,2,2);
    imagesc(AmbCall2), axis xy, colorbar
    xlabel('relative x-position (km)'), ylabel('relative y-position (km)');  
    
    % Similarity score and distance travelled
    [dist, corrScore]= crossCorrSimScores(AmbCall2, AmbCall1, ...
        deltalat_space, deltalon_space);
    
        score =[score, corrScore];
        dists=[dists, dist];
    
%     Second MR method
%     [corrScore]= crossCorrSimScoresFilt(AmbCall2, AmbCallProj, ...
%         filt_size(jj,1), filt_size(jj,2));
%     
%     score =[score, corrScore];
%     dists=[dists, dist];
    
end


figure;
plot(timegaps, score)
xlabel('delta TDOA')
ylabel('Similarity Score')
title('Two Hyd, Two Hyd')


