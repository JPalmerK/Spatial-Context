function Sim_mat= simMatMaxofProd(simStruct)
% This function creates the simulation matrix using the low
% memory approach. This should be used in most cases where not
% exploring the algorithims in depth.

% Check if the arrival table is present if not update it
if isempty(simStruct.TDOA_vals)
    disp(['Updating TDOA values'])
    UpdateTDOA(simStruct);
    
end

simStruct.titleStr ='Call Space Similarity Ideal';

% Grid X/Y space
% Get distance in meters between the lower and upper right
grid_v = vdist(min(simStruct.array_struct.latgrid),...
    min(simStruct.array_struct.longrid),...
    max(simStruct.array_struct.latgrid),...
    min(simStruct.array_struct.longrid));


% Get distance in meters between the lower left and lower right
grid_h = vdist(min(simStruct.array_struct.latgrid),...
    min(simStruct.array_struct.longrid),...
    min(simStruct.array_struct.latgrid),...
    max(simStruct.array_struct.longrid));


% Create the empty similarity matrix
arrivalArray= (simStruct.arrivalArray);

Sim_mat = zeros(size(arrivalArray,1))/0;

% Grid X/Y space
deltalat_space = (grid_v/ (length(simStruct.array_struct.latgrid)-1));
deltalon_space = (grid_h/ (length(simStruct.array_struct.longrid)-1));

% How many grid squares per second can the whale move
lat_persec = simStruct.s / deltalat_space;
lon_persec = simStruct.s / deltalon_space;


sig_tot = sqrt(simStruct.PosUncertsigma + simStruct.drift^2);
%profile on
% Step through each arrival and get it's grid probability as
% well as the projected grid probabilities for times at all
% subsiquent calls but within the maximum time cuttoff




for ii =1:length(arrivalArray)
    
    % Pick the set of call delays (number of calls)
    
    % Get the range filter
    delays = simStruct.TDOA_vals(ii, :);
    hyd_det = [simStruct.array_struct.master ...
        simStruct.array_struct.slave(simStruct.child_idx(~isnan(delays)))];
    
    
    
    % Need to send to CPU for image dialate function
    averageLklhd_space = getTruHdSpaceProd(simStruct, ii, sig_tot);
    
    % Figure out the number of time gaps within the maximum
    % allowed correlation time (time_cut)
    time_gaps = arrivalArray(ii:end, 1)-...
        arrivalArray(ii, 1);
    
    diff_times = diff(arrivalArray(ii:end, 1));
    
    % Find first big gap
    idx_end = find(diff_times>= simStruct.maxEltTime,1);
    
    if isempty(idx_end)
        idx_end = length(time_gaps);
    end
    
    time_gaps = time_gaps(1:idx_end);
    
    Lklhd_space_proj_out =  ElipsFilt(simStruct,averageLklhd_space, time_gaps,...
        grid_v,grid_h, ii);
    
    
    % If there are more than one time gap over which we need to
    % look then do the projections
    
    % Step through the time gaps, project the current likelihood surface
    % and compare it to the next likelihood surfaces
    for jj= 1:length(time_gaps)
        
        
        % Step through the time gaps/sigma values getting each
        % probability loc space and projection
        
        Lklhd_space_proj=(...
            squeeze(Lklhd_space_proj_out(:,:,jj)));
        % Get the prob. loc space space of the next call in
        % the series
        nextLklhdSpace = ...
            (getTruHdSpaceProd(simStruct,(ii+jj-1), sig_tot));
        
        
        delaysnew = simStruct.TDOA_vals(ii+jj-1, :);
        delaysnew = simStruct.array_struct.slave(simStruct.child_idx(~isnan(delays)));
        delaysnew = unique([delaysnew simStruct.array_struct.master]);
        
        % Filter the next likelihood projected space by range
        nextLklhdSpace = nextLklhdSpace .*...
            prod(simStruct.filtGrid(:,:,delaysnew),3);
        

        AScompare = prod(cat(3,Lklhd_space_proj, nextLklhdSpace),3);
        %Stic
        sim = max(AScompare(:));
        
        simValue(jj)=  sim;
    end
    
    Sim_mat(ii, ii:ii+length(simValue)-1)= simValue;
    Sim_mat(ii:ii+length(simValue)-1,ii)= simValue;
    
    simValue =[];
    disp([num2str(ii), ' of ', num2str(length(arrivalArray))])
    
end



end