function simMat1Dxcorr(obj)
% This function creates the simulation matrix using the low
% memory approach. This should be used in most cases where not
% exploring the algorithims in depth.


% Grid X/Y space
% Get distance in meters between the lower and upper right
grid_v = vdist(min(obj.array_struct.latgrid),...
    min(obj.array_struct.longrid),...
    max(obj.array_struct.latgrid),...
    min(obj.array_struct.longrid));


% Get distance in meters between the lower left and lower right
grid_h = vdist(min(obj.array_struct.latgrid),...
    min(obj.array_struct.longrid),...
    min(obj.array_struct.latgrid),...
    max(obj.array_struct.longrid));


% Create the empty similarity matrix
Sim_mat = zeros(length(obj.TDOA_vals))/0;


% Grid X/Y space
deltalat_space = grid_v/ (length(obj.array_struct.latgrid)-1);
deltalon_space = grid_h/ (length(obj.array_struct.longrid)-1);

% How many grid squares per second can the whale move
lat_persec = obj.s / deltalat_space;
lon_persec = obj.s / deltalon_space;

sig_tot = sqrt(obj.PosUncertsigma + obj.drift^2);

% Step through each arrival and get it's grid probability as
% well as the projected grid probabilities for times at all
% subsiquent calls but within the maximum time cuttoff
for ii =1:length(obj.arrivalArray)
    
    % Get the average prob loc space of the i-th call with
    % delta sigma t
    
    
    averageLklhd_space = getTruHdSpace(obj, ii, sig_tot);
    
    % Figure out the number of time gaps within the maximum
    % allowed correlation time (time_cut)
    time_gaps = obj.arrivalArray(ii:end, 1)-...
        obj.arrivalArray(ii, 1);
    time_gaps = time_gaps(time_gaps<obj.time_cut);
    
    % If there are more than one time gap over which we need to
    % look then do the projections
    
    if length(time_gaps)>1
        
        % Calculate the sigma values based on EM Nosal
        % suggestion
        % sigma (error) valuse from the normal distribution
        
        
        % Step through the time gaps/sigma values getting each
        % probability loc space and projection
        for jj= 1:length(time_gaps)
            
            
            nextLklhdSpace = getTruHdSpace(obj, (ii+jj-1), sig_tot);
            
            % Get the comparison value of the projected prob
            % loc space and the next call in the squence
            aa = sum(averageLklhd_space,1);
            bb = sum(averageLklhd_space,2);
            
            
            aa2 = sum(nextLklhdSpace,1);
            bb2 = sum(nextLklhdSpace,2);
            
            cor1 = corrcoef(aa,aa2);
            cor2 = corrcoef(bb,bb2);
            
            simValue = (cor1(1,2) + cor2(1,2))/2;
            
            
            % Populate the simulation matrix
            Sim_mat(ii, ii+jj-1) = simValue;
            Sim_mat(ii+jj-1 ,ii) = simValue;
            
        end
        
    end
    %disp([num2str(ii), ' of ',...
    %    num2str(length(obj.arrivalArray))])
end
obj.Sim_mat= Sim_mat;



end