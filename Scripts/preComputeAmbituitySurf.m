function propAmbSurf = preComputeAmbituitySurf(simStruct, fname, MaxTimeMin )


% Maximum time beyond which to not project calls
if nargin < 2
    fname = 'ReginaAllDet.mat';
end


% Maximum time beyond which to not project calls
if nargin < 3
    MaxTimeMin = 15;
end
parent = simStruct.array_struct.master;
MaxTime = MaxTimeMin*60;

sig_tot = sqrt(simStruct.PosUncertsigma + simStruct.drift^2);

% Pull out the arrival array for readability
arrivalArray= (simStruct.arrivalArray);

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
close all
for ii=1:size(arrivalArray,1)
    
    
    delays = simStruct.TDOA_vals(ii, :);
    
    if all (isnan(delays))
        disp('something broke')
    else
        
        
        hyd_det = [simStruct.array_struct.master ...
            simStruct.array_struct.slave(simStruct.child_idx(~isnan(delays)))];
        filt = prod(simStruct.filtGrid(:,:,hyd_det),3);
        
        % Need to send to CPU for image dialate function
        averageLklhd_space = getTruHdSpaceProd(simStruct, ii, sig_tot);
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

parfor ii=1:size(arrivalArray,1)
    
    propAmbSurf = struct();
    idxs = simStruct.arrivalTable.dex(ii:end);
    
    
    fname = strcat([simStruct.propAmbLoc, '\Parent_',num2str(parent),...
        '_Dex_', num2str(idxs(1)), '.mat']);
    if ~isfile(fname)
        
        disp(['computing', num2str(ii), ' of ' num2str(size(arrivalArray,1))])
        delays = simStruct.TDOA_vals(ii, :);
        hyd_det = [simStruct.array_struct.master ...
            simStruct.array_struct.slave(simStruct.child_idx(~isnan(delays)))];
        
        %     Need to send to CPU for image dialate function
        averageLklhd_space = getTruHdSpaceProd(simStruct, ii, sig_tot);
        
        %     Figure out the number of time gaps within the maximum
        %     allowed correlation time (time_cut)
        time_gaps = arrivalArray(ii:end, 1)-...
            arrivalArray(ii, 1);
        
        diff_times = [0; diff(arrivalArray(ii:end, 1))];
        
        %     Don't bother looking for calls beyond 15 minutes
        
        idxs = idxs(time_gaps<MaxTime);
        time_gaps = time_gaps(time_gaps<MaxTime);
        
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