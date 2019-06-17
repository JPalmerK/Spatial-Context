
function [chain, Cluster_id] = MultiplHydrophonClusterDelays(array_struct,...
    localize_struct,array_idx,fs, NChidlHyd, time_cut,...
    cor_thresh, LSQ_plot_thres, hydrophone_struct)
% Function to group detections into clusters


% % Inputs:
%     array_struct - array structure from GPL
%     localize_struct - structure of localized detections from GPL
%     array_idx - which sub array (for DCLDE use 1, for PMRS use 19)
%     fs - sample frequency
%     NChildHyd - number of pairs of hydrophones used to make the final corr plot
%     time_cut - maximum time (min) that calls can be correlated over
%     corr_thres - minimum correlation coefficient used to link calls
%     LSQ_plot_thres - for plotting only, minimum LSQ threshold a call must have to be plotted
%     hydrophone_struct - optional parameter, plot hydrophone locations if
%     available
% Returns:
%    chain - chain structure created by make_Ladder_linkages
%    Cluster_id - array length number of detections containing the group
%    ide for each call


%% Look at the matches on each parent/child pair

% Get the parent hydrophone
ParentHyd = array_struct(array_idx).master;

% Get the children hydrophones
ChildHyds = array_struct(array_idx).slave;



%Create the correlation plots for each parent/child pair
Corr_coef_map = zeros(NChidlHyd,...
    length(localize_struct.hyd(ParentHyd).delays),...
    length(localize_struct.hyd(ParentHyd).delays));


% max swim speed
s =12;
c=1500;
% For each pair of hydrophones calculate the potential TDOA space for all
% the calls



for ii=1:1
    
    % Get the delays for all the calls
    for kk = 1:size(localize_struct.hyd(ParentHyd).delays,1)
        
       
        % potential toa space
        deltatdoa = (localize_struct.hyd(ParentHyd).delays(:,ii) -...
            localize_struct.hyd(ParentHyd).delays(kk,ii)).^2;
        
        % Time between call k and other calls in seconds (delays in seconds)
        deltat = (localize_struct.hyd(ParentHyd).rtimes - ...
                   localize_struct.hyd(ParentHyd).rtimes(kk))/fs;
        
        %sigma = (s/c)*min(abs(deltat), 5*60);
        %sigma = (s/c)*abs(deltat);
        sigma = (deltat * (2*s)/c);
        
        
        
        % Grow the LSQ space
        LSQ_space =  normpdf(deltatdoa, 0, abs(1+sigma(kk)));    

        % Normalize it
        LSQ_space = (LSQ_space - min(min(LSQ_space))) /...
            ( max(max(LSQ_space)) - min(min(LSQ_space)));
            
                
        % Squash it
        Corr_coef_map(ii,kk,:) = LSQ_space;
        % Now we have to sum of the TDOA space for each call
        
    end
    
    
    
end
imagesc(squeeze(Corr_coef_map(1,:,:)))
% Make correlation coefficient maps for the reference pair and all the
% other pairs
nrows = floor(sqrt(NChidlHyd));
ncols = ceil(NChidlHyd/nrows);

% If using more than one hydrophone pair, get the mean correlation across hydrophone
if NChidlHyd>1
    Corr_coef_map_out = squeeze(mean(Corr_coef_map));
else
    Corr_coef_map_out = squeeze(Corr_coef_map);
end

figure
imagesc(Corr_coef_map_out)

%% Make ladder structures and plot

% correlate out to 15 min (dive time of a whale)

% time cut in seconds
time_cut = time_cut*60;


% Set up initial times
if length(size(localize_struct.hyd(array_struct(array_idx).master).coordinates))>2
    x=localize_struct.hyd(array_struct(array_idx).master).coordinates(2,2,:);
    y=localize_struct.hyd(array_struct(array_idx).master).coordinates(2,1,:);
    score=localize_struct.hyd(array_struct(array_idx).master).score(2,:);
    
else % it's a simulation so don't bother for now
    x=localize_struct.hyd(array_struct(array_idx).master).coordinates(:,2);
    y=localize_struct.hyd(array_struct(array_idx).master).coordinates(:,1);
    score = localize_struct.hyd(array_struct(array_idx).master).score;
end

times=localize_struct.hyd(array_struct(array_idx).master).rtimes;

% calltype=localize_struct.hyd(array_struct(j).master).calltype;

x=squeeze(x);
y=squeeze(y);
score=squeeze(score);


% Creat the chain links
chain = make_Ladder_linkeges(Corr_coef_map_out,...
    localize_struct, ParentHyd, fs, cor_thresh, time_cut);

%%  Make big datafram ewith all the relevent information
% Call id's for each call
Cluster_id = ones(size(x));

for ii=1:length(chain);
    Cluster_id(chain(ii).index) = ii;
    
end
%% Make Plots

% If this is a simulation color plots by group, not arrival time
if isfield(localize_struct.hyd(ParentHyd), 'agent_ids')
    
    
    title1 =[ 'Data Simulation Results:' newline...
        'Predicted Clusters Using Correlation Between ', ...
        num2str(NChidlHyd),' Hydrophone Pairs'];
    title2 = 'True Clusters';
    
    % Number of unique colors needed
    ncolors = max([length(unique(Cluster_id))...
        length(unique(localize_struct.hyd(ParentHyd).agent_ids))]+1);
    ColorVals = lines(ncolors);
    
    
    % Color scales for simulated clusters
    color_var_drift = ColorVals(Cluster_id,:);
    colorscale2 = ColorVals(localize_struct.hyd(ParentHyd).agent_ids,:);
    
    % Color bar titles
    cbar1label = 'Predicted Clusters';
    cbar2label = 'True Clusters';
    
    x_trim_drift= x;
    y_trim_drift =y;
    colormap(lines(ncolors))
    
    txt_str2 =[num2str(length(unique(localize_struct.hyd(ParentHyd)...
        .agent_ids))), ' Agents'];
    txt_str1 =[num2str(length(unique(Cluster_id))), ' Groups Predicted'];
    x_txt = min(x)+ range(x)/2;
    y_txt = min(y); 
    


else % real data, clusters unknown plot by arrival times
    % trim the localizations if min LSQ provided
    if exist('LSQ_plot_thres','var') == 1
        [k1 k2]= find(score <LSQ_plot_thres);
    else
        k2 = 1:length(x);
    end
    
    title1 =[ 'Localized Data Clusterd Using Ladder Structure and Mean'...
        newline 'Correlation Between ', num2str(NChidlHyd),' Hydrophone Pairs'];
    title2 = 'Localized Data Plotted Against Time of Detection';
    
    color_var_drift= zeros(size(x));
    ColorVals = jet(length(k2));
    colorscale2 = (localize_struct.hyd(ParentHyd).rtimes(k2)...
        - min(localize_struct.hyd(ParentHyd).rtimes(k2)))/fs;
    

    
    % Create color variable
    cluster_n = 1; % initial cluster color (other than 0)
    
    % Color all localizations in the cluster the same
    for ii=1:length(chain)
        
        if chain(ii).n >1
            color_var_drift([ii, chain(ii).index]) = (log10(chain(ii).n))+.001;
            cluster_n =cluster_n+1;
        end
    end
    

    
    x_trim_drift=x(k2);
    y_trim_drift=y(k2);
    color_var_drift = color_var_drift(k2);
    
    
    cbar1label = 'log10(Number of Detections in Cluster)';
    cbar2label = 'Elapsed Time Since First Detection (sec)';
    
    % text labels, none for unknown data
    txt_str2 =[''];
    x_txt = max(x)-.1
    y_txt = min(y) +.01
    
    
    txt_str1 =[num2str(length(unique(Cluster_id))), ' Groups Predicted'];
    y_txt = min(y);  
    
end

figure
fig = get(groot,'CurrentFigure');

%Plot data
fig
subplot(2,1, 1)
scatter(x_trim_drift, y_trim_drift, 20, color_var_drift, 'filled');
title(title1)
ylabel('Lat'); xlabel('Lon')
if ~isempty(LSQ_plot_thres) 
    colormap((jet))
else
    colormap(lines(ncolors))
end
text(x_txt, y_txt, txt_str1)
h = colorbar;
ylabel(h, cbar1label)
hold off

fig
subplot(2,1, 2)
scatter(x_trim_drift, y_trim_drift, 20, colorscale2, 'filled');
title(title2)
text(x_txt, y_txt, txt_str2)
ylabel('Lat'); xlabel('Lon')
h = colorbar;
ylabel(h, cbar2label)



% If they hydrophone structure is present (e.g. not PMRF) add the
% hydrophones
if ~isempty(hydrophone_struct)


    hyd_arr = struct2cell(hydrophone_struct);
    hyd_arr =vertcat(hyd_arr{2,:,:});

    subplot(2,1, 1)
    hold on
    scatter(hyd_arr(:,2), hyd_arr(:,1), 100, 'k', 'filled', '^');


    subplot(2,1, 2)
    hold on
    scatter(hyd_arr(:,2), hyd_arr(:,1), 100, 'k', 'filled', '^');



end










%%
% %drifts = [-10 -2 -.5 0 .5 2 10];
% drifts = linspace(-1,1,11);
%
%
% if ~exist('RefID','var')
%     % third parameter does not exist, so default it to something
%     RefHydID = array_struct(1).slave(1)
% end
%
% SlaveHydID = array_struct(1).slave
% slave_idx= 1:length(SlaveHydID);
%
%
% Px = zeros(length(array_struct(1).latgrid)*...
%     length(array_struct(1).longrid), 1);
% % Number of pairs of hydrophones to look at to include reduce
% % localization error
%
%
% % Get the delays for all the calls
% for kk = 1:size(localize_struct.hyd(19).delays,1)
%
%     % Pick a set of call delays (number of calls)
%     delays = localize_struct.hyd(19).delays(kk,:);
%
%     % Get the tdoa space
%     toa_space = (cell2mat(array_struct(1).toa_diff(slave_idx(1)+1)));
%
%     % potential toa space
%     LSQ_VAL = abs(toa_space - delays(slave_idx(1)));
%
%     %LSQ_space = LSQ_VAL;
%     LSQ_space = 1./(1+exp(LSQ_VAL));
%     %figure; imagesc(LSQ_space)
%
%     % Squash it
%     Px(:,kk) = reshape(LSQ_space.',[] ,1);
%     % Now we have to sum of the TDOA space for each call
%
% end
%
% %figure(1)
% Corr_coef_map_ref = corrcoef(Px);
% % imagesc(Corr_coef_map_ref)
% %title('Refrence Correlation Coefficent Map H1 and H2')
%
%
% % Normalize all values
% Px = ((Px - min(Px)) ./ ( max(Px) - min(Px)));
%
%
% %
% % Px1 = zeros(length(array_struct(1).latgrid)*...
% %     length(array_struct(1).longrid), 1);
% % Corr_coef_map = zeros(length(drifts),...
% %     length(localize_struct.hyd(19).delays),...
% %     length(localize_struct.hyd(19).delays));
% %
% %
% % for zz=1:length(drifts)
% %     Get the delays for all the calls
% %     for kk = 1:size(localize_struct.hyd(19).delays,1)
% %
% %         Pick a set of call delays (number of calls)
% %         delays = localize_struct.hyd(19).delays(kk,:);
% %
% %
% %         Get the tdoa space
% %         toa_space = (cell2mat(array_struct(1).toa_diff(3)));
% %
% %         potential toa space
% %         LSQ_VAL = abs(toa_space - (delays(slave_idx(2)))+drifts(zz));
% %
% %
% %         LSQ_VAL = LSQ_VAL - mean(mean(LSQ_VAL));
% %         LSQ_space = 1./(1+exp(LSQ_VAL));
% %         figure; imagesc(LSQ_space)
% %
% %         Squash it
% %         Px1(:,kk) = reshape(LSQ_space.',[] ,1);
% %         Now we have to sum of the TDOA space for each call
% %
% %     end
% %
% %
% %     Create the correlation heat map
% %     Corr_coef_map(zz,:,:) = corrcoef(Px);
% %     disp(num2str(zz))
% %     figure(2)
% %     subplot(3,3,zz)
% %     imagesc(squeeze(Corr_coef_map(zz,:,:)))
% %     title(['H1 H3 Correlation with ' num2str(drifts(zz)) ' sec Drift'])
% % end
%
%
% %% Try fuzzy logic approach
%
%
% % Look at the standard deviation of the values greater than 0, low std
% % suggests more focused lense
%
% close all
%
% %   Refine the space to only calls with high correlation scores on both
% %  hydrophones
%
% % Index of the scores greater than .8 on the first hydrophone pairs
% idx1 = find(localize_struct.hyd(19).cross_score(:,1)>0.7);
% idx2 = find(localize_struct.hyd(19).cross_score(:,1)>0.7);
%
% idx_keep =find(all(~isnan(localize_struct.hyd(19).cross_score), 2));
%
% outscore = zeros(size(drifts));
% area_under_curve_vals = zeros(length(idx_keep), length(drifts));
%
% for zz=1:length(drifts)
%
%     % Get the delays for all the calls
%     for kk = 1:size(localize_struct.hyd(19).delays,1)
%
%         % Pick a set of call delays (number of calls)
%         delays = localize_struct.hyd(19).delays(kk,:);
%
%
%         % Get the tdoa space
%         toa_space = (cell2mat(array_struct(1).toa_diff(3)));
%
%         % potential toa space
%         LSQ_VAL = abs(toa_space - (delays(slave_idx(2)))+drifts(zz));
%
%         %
%         %         LSQ_VAL = LSQ_VAL - mean(mean(LSQ_VAL));
%         LSQ_space = 1./(1+exp(LSQ_VAL));
%         %figure; imagesc(LSQ_space)
%
%         % Squash it
%         Px1(:,kk) = reshape(LSQ_space.',[] ,1);
%         % Now we have to sum of the TDOA space for each call
%
%     end
%
%
%     % Normalize all values
%     Px1 = ((Px1 - min(Px1)) ./ ( max(Px1) - min(Px1)));
%
%
%     % Fuzzy logic and the two spaces
%     mm = Px+Px1;
%
%     figure(4)
%
%     subplot(length(drifts), 3, (zz*3)-2)
%     imagesc(reshape(Px(:,5), 301,216))
%
%     subplot(length(drifts), 3, (zz*3)-1)
%     imagesc(reshape(Px1(:,5), 301,216))
%
%
%     subplot(length(drifts), 3, zz*3)
%     imagesc(reshape(mm(:,5), 301,216))
%     title([num2str(drifts(zz)) ' s Offset ' num2str(sum(mm(:,5))/ length(mm)) ])
%
%
%     figure(5)
%
%
%     subplot(length(drifts), 3, (zz*3)-2)
%     imagesc(reshape(Px(:,1), 301,216))
%
%     subplot(length(drifts), 3, (zz*3)-1)
%     imagesc(reshape(Px1(:,1), 301,216))
%
%
%     subplot(length(drifts), 3, zz*3)
%     imagesc(reshape(mm(:,1), 301,216))
%     title([num2str(drifts(zz)) ' s Offset ' num2str(sum(mm(:,1))/ length(mm)) ])
%
%
%
%
%     % Get the standard deviation of all points where the value is greater
%     %     % than 0
%     %
%     %     figure(5)
%     %     subplot(5,1,zz)
%     %
%     %     hist(reshape(mm(mm>0),[],1))
%     %
%     %     title(['Drift Shifted By ' num2str(drifts(zz)) ' seconds'])
%     %     val(jj)= std(aa(2:end)-aa1);
%     %     text(-.25,300, ['St.Dev ' num2str(val(jj))])
%     %     ylim([0,400])
%     %     xlim([-.3 .3])
%
%
%
%
%     area_under_curve_vals(:,zz) = sum(mm(:,idx_keep))/ length(mm(:,idx_keep));
%     % Total number of squares in the 95th percentile
%     percentile_vals = mean(sum(mm(:,idx_keep))/ length(mm(:,idx_keep)));
%
%     %     % total number of squares in that percentile
%     %     aa = (mm>1);
%     %
%     %     % Ratio of non zeros
%     %     aa = percentile_vals ./ sum(aa);
%     %
%     outscore(zz) = mean(percentile_vals);
%
% end
%
%
%
%
%
% p = polyfit(drifts,outscore, 4)
% x1 = linspace(min(drifts), max(drifts),20);
% y1 = polyval(p,x1);
% [~, aa] = max(y1)
% figure
% plot(x1, y1)
% hold on
% scatter(drifts, outscore)
% text(mean(x1), mean(y1), ['min value at ' num2str(x1(aa)) ' s offset'])
% xlabel('Simulated Clock Drift')
%
%
% % determine if any of the max values are on the edge
%
% for zz=1:length(drifts)
% end
%
%
%
%
% figure(3)
% for zz=1:length(drifts)
%
%     subplot(4,4,zz)
%     hist(area_under_curve_vals(:,zz),20)
%     title([num2str(drifts(zz)) ' s Offset '])
%     xlim([0 .8])
%     ylim([0 100])
%
%
% end
%
% subplot(4,4,16)
% plot(drifts,  quantile(area_under_curve_vals, .5))
% hold on
% line([0 0], [.41 .43]);
% hold off
%
%
% figure
%
% %
% %
% % p = polyfit(drifts, quantile(area_under_curve_vals, .5), 2)
% % x1 = linspace(min(drifts), max(drifts),20);
% % y1 = polyval(p,x1);
% % [~, aa] = max(y1)
% % figure
% % plot(x1, y1)
% % hold on
% % scatter(drifts, quantile(area_under_curve_vals, .5))
% % text(mean(x1), mean(y1), ['min value at ' num2str(x1(aa)) ' s offset'])
% % xlabel('Simulated Clock Drift')
%
%
%
% %%
%
% Fs = 6000;  % magic number
% aa =[]; bb=aa;
% val = zeros(length(drifts),1);
%
% % Corr_coef_map contains the TDOA space for a pair of hydrophones
% % Loop through
%
%
% for jj = 1:length(drifts)
%
%     meanimage = abs((Corr_coef_map_ref+squeeze(Corr_coef_map(jj,:,:))));
%     meanimage = 1./(1+exp(-meanimage));
%     for ii=1:2
%
%
%         aa1 = reshape(meanimage(ii,ii+1:end),[],1)';
%
%         time_breaks = (localize_struct.hyd(19).rtimes(ii+1:end)/Fs)';
%
%         bb1=  time_breaks; %(time_breaks - time_breaks(1))';
%
%
%         cc = repmat(ii, 1, length(bb1));
%         figure(3)
%         subplot(7,1,jj)
%         scatter(bb1, aa1, 20, cc, 'filled')
%         ylabel('Correlation Values')
%         title(['Drift Shifted By ' num2str(drifts(jj)) ' seconds'])
%         hold on
%
%         if jj ==5
%             xlabel('Time Difference')
%         end
%
%         if ii==2
%             figure(5)
%             subplot(7,1,jj)
%
%             hist(aa(2:end)-aa1)
%             title(['Drift Shifted By ' num2str(drifts(jj)) ' seconds'])
%             val(jj)= std(aa(2:end)-aa1);
%             text(-.25,300, ['St.Dev ' num2str(val(jj))])
%             ylim([0,400])
%             xlim([-.3 .3])
%         end
%
%
%         aa = [aa1];
%         bb = [bb1];
%     end
%
%     figure(4)
%     subplot(3,3,jj)
%     imagesc((meanimage), [.6 .85])
%     title(['Drift Shifted By ' num2str(drifts(jj)) ' seconds'])
%
%
% end
%
%
%
%
%
%
%
%
%
%
% %%
%
% % Create entropy structure
% entropy_measure = zeros(length(drifts), 7);
%
% [r,c] = size(Corr_coef_map_ref);
% n = sum(sum(Corr_coef_map_ref))/2;
% p = Corr_coef_map_ref/n;  % probability estimator
%
%
% for ii=1:length(entropy_measure)
%
%     totalP=0;
%     meanimage = abs((Corr_coef_map_ref+squeeze(Corr_coef_map(ii,:,:))));
%     meanimage = 1./(1+exp(-meanimage));
%
%     % Get the time breaks between the calls
%     time_breaks = diff((localize_struct.hyd(array_id).rtimes/6000))';
%
%     breaks_idx = find(time_breaks>900);
%     %
%     %     % Convert to probabilities
%     %     meanimage = meanimage / sum(sum(meanimage));
%
%     %
%     %
%     %     figure(1)
%     %     subplot(3,3,ii)
%     %imagesc(meanimage)
%     aa = hist(diff(meanimage(:,1)), [-.5:.01:0.5]);
%
%     hist(diff(meanimage(:,1)),[-.5:.01:0.5] );
%     %     xlabel(['Seconds added ' num2str(drifts(ii))])
%     %     temp = [];
%     %     temp_tot =0;
%     %[entropy_measure(ii,1) entropy_measure(ii,2)]= entropy(meanimage(find(triu(meanimage,1))));
%     entropy_measure(ii,3) = std(abs(diff(meanimage(:,1))));
%     entropy_measure(ii,4) = iqr(aa);
%     entropy_measure(ii, 5) = prctile((abs(diff(meanimage(:,1)))), 5);
%     %     % calculate an entropy measure - rough area round the center of the mat
%     %     for jj=1:size(Corr_coef_map_ref,2);
%     %
%     %
%     %
%     %         if ~ismember(jj, breaks_idx)
%     %             % get the end of the block
%     %             block_end = breaks_idx(find(breaks_idx>jj,1));
%     %             % Get the miniblock
%     %             calls_of_interest = meanimage(jj+1:block_end,jj);
%     %             % NOrmalize it
%     %             p2 = calls_of_interest / (sum(calls_of_interest));
%     %
%     %             % Calculate the entropy of the row
%     %             totalP = totalP + (-sum(p2.*log2(p2)));
%     %
%     %             temp = [temp calls_of_interest'];
%     %
%     %           end
%     %
%     %     length(temp)
%     %     end
%     %
%     %
%     %     entropy_measure(ii) = totalP;
%     %     entropy_measure1(ii) = (-sum(temp.*log2(temp)));
%     %
%
%
%
%     % Try something else, plot the correlation vs lag times between calls
%     % and fit a gaussian mixture model
%     X = [reshape(meanimage,[],1), log10(repmat(time_breaks, length(meanimage)))];
%     bic =[];
%
%
%
%
%     for jj = 10:20
%
%         GMModel = fitgmdist(X,jj,'RegularizationValue',.0001);
%         bic =[bic GMModel.BIC];
%
%         %         figure(4)
%         %         h = gscatter(X(:,1),X(:,2));
%         %         hold on
%         %         ezcontour(@(x1,x2)pdf(GMModel,[x1 x2]),get(gca,{'XLim','YLim'}))
%         %         title('{\bf Scatter Plot and Fitted Gaussian Mixture Contours}')
%         %         legend(h,'Model 0','Model1')
%         %         hold off
%
%     end
%     [cc,xx] = min(bic);
%     entropy_measure(ii,6) = cc;
%
%
% end
% entropy_measure(:,7) = drifts;
% figure(2)
% plot(drifts, abs(entropy_measure(:,6)), '-.k*');
% xlabel('Standard Deviation of Diff Histogram')
% ylabel('Image Entropy')


