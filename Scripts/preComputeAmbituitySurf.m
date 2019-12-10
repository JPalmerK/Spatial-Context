function propAmbSurf = preComputeAmbituitySurf(simStruct, fname, MaxTimeMin )


% Maximum time beyond which to not project calls
if nargin < 3
    MaxTimeMin = 15;
end


parent = simStruct.array_struct.master;
MaxTime = MaxTimeMin*60;

% Position uncertainty 
sig_tot = sqrt(simStruct.PosUncertsigma + simStruct.drift^2);

% Pull out the arrival array for readability
arrivalArray= (simStruct.arrivalArray);

% Grid granularity
grid_v = vdist(min(simStruct.array_struct.latgrid),...
    min(simStruct.array_struct.longrid),...
    max(simStruct.array_struct.latgrid),...
    min(simStruct.array_struct.longrid));


% Get distance in meters between the lower left and lower right
grid_h = vdist(min(simStruct.array_struct.latgrid),...
    min(simStruct.array_struct.longrid),...
    min(simStruct.array_struct.latgrid),...
    max(simStruct.array_struct.longrid));


% Make an initial array containing all the ambiguity surfaces propagated to
% Parent_5_Dex_4934
initAmbSurfs = zeros([length(simStruct.array_struct.latgrid),...
    length(simStruct.array_struct.longrid), size(arrivalArray,1)]);


% Create ambiguity surfaces for all calls at the time they arrived
% (projected time = 0)
close all
for ii=1:size(arrivalArray,1)
    
    % TDOA values for each call
    delays = simStruct.TDOA_vals(ii, :);
    
    if all (isnan(delays))
        disp('something broke')
    else
        
        % Hydrophone numbers where the call was detected 
        hyd_det = [simStruct.array_struct.master ...
            simStruct.array_struct.slave(simStruct.child_idx(~isnan(delays)))];
        
        % Pull out the filter grid (ie detection function) for each of the
        % hydrophones where the call was detected. Take the product to
        % identify the football shaped area from where the call could have
        % originated 
        filt = prod(simStruct.filtGrid(:,:,hyd_det),3);
        
        % Create the ambiguity surface with the sim structure, the index of
        % the call and the total uncertainty 
        averageLklhd_space = getTruHdSpaceProd(simStruct, ii, sig_tot);
        
        % Multiply ambiguity surface by the deteciton function surface to
        % remove distance areas where calls could have origniated (only
        % really usefull when calls are detected by 2 hydrophones only
        averageLklhd_space= (averageLklhd_space).*filt;
        initAmbSurfs(:,:,ii)=averageLklhd_space;
    end
    
end

fname = strcat([simStruct.propAmbLoc,'\InitialAmbSurfParent_'...
    num2str(parent)]);

AmbSurfs = struct();
AmbSurfs.ProjTimeIdxs = simStruct.arrivalTable.dex;
AmbSurfs.AmpSurfs = initAmbSurfs;

save(fname,'AmbSurfs', '-v7.3')

%
% Stopped at 3375, continue later Parent_5_Dex_4934 missing

% Loop through each call and create the ambiguity surfaces and projected
% ambituity surfaces. Save these along with the call/detection index (dex
% from localize structure) and the times of the subsiquent calls for which
% the call was projected

parfor ii=1:size(arrivalArray,1)
    
    propAmbSurf = struct();
    idxs = simStruct.arrivalTable.dex(ii:end);
    
    % File name where to store the projected ambiguity surfaces
    fname = strcat([simStruct.propAmbLoc, '\Parent_',num2str(parent),...
        '_Dex_', num2str(idxs(1)), '.mat']);
    
    % Check if it already exists, this line needs to be removed if any
    % parameters are changed
    if ~isfile(fname)
        
        disp(['computing', num2str(ii), ' of ' num2str(size(arrivalArray,1))])
        delays = simStruct.TDOA_vals(ii, :);
        hyd_det = [simStruct.array_struct.master ...
            simStruct.array_struct.slave(simStruct.child_idx(~isnan(delays)))];
        
        %     Need to send to CPU for image dialate function
        averageLklhd_space = getTruHdSpaceProd(simStruct, ii, sig_tot);
        
        % Find the time difference between the current call (ii) and all
        % the subsiquent calls in the detection file
        time_gaps = arrivalArray(ii:end, 1)-...
            arrivalArray(ii, 1);
        
        % Determine the time differences between the call and all
        % subsiquent calls (look to start a new acoustic encounter(
        diff_times = [0; diff(arrivalArray(ii:end, 1))];
        
        %     Don't bother looking for calls beyond 15 minutes
        idxs = idxs(time_gaps<MaxTime);
        
        time_gaps = time_gaps(time_gaps<MaxTime);
        
        % Create the projected ambiguity surfaces for all arrival times not
        % knocked out in the above code
        Lklhd_space_proj_out =  ElipsFilt(simStruct,averageLklhd_space,...
            time_gaps, grid_v,grid_h, ii);
        
        Lklhd_space_proj_out(Lklhd_space_proj_out< .001) = 0;
        
        sparseOut = struct();
        for jj=1:size(Lklhd_space_proj_out, 3)
            sparseOut(jj).ambSurf =sparse(Lklhd_space_proj_out(:,:,jj));
        end
        
        %  Populate the projection structure using the indexes from the GPL
        % detections
        propAmbSurf.AmpSurfs =  sparseOut;
        propAmbSurf.ProjTimeIdxs = idxs;
        propAmbSurf.deltaSec =time_gaps;
        parsave(fname, propAmbSurf)
        disp('File Written')
        
    end
end





end