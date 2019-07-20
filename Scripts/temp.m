% Run experiments
close all; clear all
clear classes; clc

whereAmI = loadSimspace();
dclde_2013_meta = xlsread(whereAmI{1});
load(whereAmI{2});

sigma = 0.0004^2 +.1^2 + .3^2;
grid= cell2mat( array_struct.toa_diff(2));
grid2 =  cell2mat( array_struct.toa_diff(3));

%% Example 1
% two well localized calls

tdoaObs = [0 3];
gridAdj1 = normpdf(grid-tdoaObs(1), 0, 2*sigma)./normpdf(0,0,2*sigma);
gridAdj2 = normpdf(grid2-tdoaObs(2), 0, 2*sigma)./normpdf(0,0,2*sigma);
gridCall1 = cat(3,gridAdj1, gridAdj2);
gridCall1 = min(gridCall1,[],3);
imagesc(gridCall1)

% second call, shifted a bit
tdoaObs =[0 2];
gridAdj1 = normpdf(grid-tdoaObs(1), 0, 2*sigma)./normpdf(0,0,2*sigma);
gridAdj2 = normpdf(grid2-tdoaObs(2), 0, 2*sigma)./normpdf(0,0,2*sigma);
gridCall2 = cat(3,gridAdj1, gridAdj2);
gridCall2 = min(gridCall2,[],3);
figure; imagesc(gridCall2)


aa =normxcorr2(gridCall2,gridCall1)
[ypeak, xpeak] = find(aa==max(aa(:)));

% correlation value 
corr_val =aa(ypeak, xpeak)
[call1y, call1x] = find(gridCall1==max(gridCall1(:)));
[call2y, call2x] = find(gridCall2==max(gridCall2(:)));

% Distance travelled
dist = sqrt((call1y-call2y)^2 +(call1x-call2x)^2)
% Speed
speed = dist/60; % assume duration between calls was 60 sec

% speed normalized likelihood
pd = makedist('Normal',4,5)
t = truncate(pd,0,12);
plot([-10:.1:10], pdf(t, [-10:.1:10])/pdf(t,4))
pdf(t, 10)/pdf(t,4)

