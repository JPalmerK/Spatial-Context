function Sim_mat= simMatMaxofProdPreComputed(simStruct)
% Create the similarity matrix for the ambiguity spaces using pre-computed
% propagatied ambiguity


% % This section is currently commented/uncommented between runs
% % Check if pre-computed ambiguity surface exists, if it's not there make
% % one
% if isfield(simStruct, 'propAmbLoc')
%     
%     try
%         load(simStruct.propAmbLoc);
%     catch
%         disp('No projected ambiguity surface grids available, computing')
%        preComputeAmbituitySurf(simStruct);
%     end
% else
%     propAmbSurf = preComputeAmbituitySurf(simStruct);
%     
% end

% Create the empty similarity matrix
arrivalArray= (simStruct.arrivalArray);
Sim_mat = zeros(size(arrivalArray,1))/0;
parent = simStruct.array_struct.master;

% Location where ambiguity surfaces are stored
floc = simStruct.propAmbLoc;

gplId = simStruct.arrivalTable.dex;
%Confirm all ambiguity surfaces are present
for ii =1:size(arrivalArray,1)
    % GPL Ids of the call

    % File name
    fname = strcat([floc, '\Parent_',num2str(parent),...
        '_Dex_', num2str(gplId(ii)), '.mat']);
    
    if ~isfile(fname)
        disp(['GPL detection ' num2str(gplId(ii)), ' missing'])
    end
end



% Load the initial ambiguity surface structure
load(strcat([floc, '\InitialAmbSurfParent_', num2str(parent),'.mat']))


% simMatTemp = Sim_mat;
% 
% for ii =1:size(arrivalArray,1)
%     gplId = simStruct.arrivalTable.dex(ii);
%     idx1 = find(AmbSurfs.ProjTimeIdxs==gplId);
%     currGrid = AmbSurfs.AmpSurfs(:,:,idx1);
%     
%     for jj = ii:size(arrivalArray,1)
%         gplId2 = simStruct.arrivalTable.dex(jj);
%         idx2 = find(AmbSurfs.ProjTimeIdxs==gplId2);
%         nextGrid = AmbSurfs.AmpSurfs(:,:,idx2);
%         val = immse(currGrid, nextGrid);
%         simMatTemp(ii,jj) = val;
%         simMatTemp(jj,ii)= val;
%         
%     end
% end
%         





disp('blarg')




% Create the similarity matrix from the projected ambiguity surfaces
for ii =1:size(arrivalArray,1)

    % GPL Ids of the call
    gplId = simStruct.arrivalTable.dex(ii);
    
    % File name
    fname = strcat([floc, '\Parent_',num2str(parent),...
        '_Dex_', num2str(gplId(1))]);
    
    % Load the ambiguity surfaces
    load(fname)
    
    % Rename this so we can load the other file later on
    projAmbSurfCurrent = x;
    
    % get the ID's for the times to proejct the calls
    projIdxs = projAmbSurfCurrent.ProjTimeIdxs;
    
    % Get the call times for the projected calls
    projIdxstrim = intersect(projIdxs,...
        simStruct.arrivalTable.dex(ii:end));
    
    simValue = zeros(length(projIdxstrim),1);
    AScompare_out = zeros([size(simStruct.filtGrid,1)...
        size(simStruct.filtGrid,2)...
        length(projIdxstrim)]);
    
    aa =full(struct2array(projAmbSurfCurrent.AmpSurfs(1)));
    for jj= 1:length(projIdxstrim)
        

        % Index value of the next ambiguity surface(not projected)
        nextDex =projIdxstrim(jj);
        
        % Index in the non-time expanded ambiguity surfaces
        idx = find(AmbSurfs.ProjTimeIdxs ==nextDex);
        
        if isempty(idx)
            disp(['Missing Null Ambiguity surf ', num2str(nextDex)]);
        end
        nextLklhdSpace = AmbSurfs.AmpSurfs(:,:,idx);
        
        
        % Get the current ambiguity surface (projected)
        projIdx = find(projAmbSurfCurrent.ProjTimeIdxs == nextDex);
        Lklhd_space_proj=full(struct2array(...
            projAmbSurfCurrent.AmpSurfs(projIdx)));
        
        AScompare = prod(cat(3,Lklhd_space_proj, nextLklhdSpace),3);
        AScompare_out(:,:,jj) =AScompare;
%         
              
%         sim = ssim(aa,nextLklhdSpace);
%         
        
        x1dist = sum(aa,1);
        y1dist = sum(aa,2);
        
        
        x2dist = sum(nextLklhdSpace,1);
        y2dist = sum(nextLklhdSpace,2);
        
        
        sim = 1-immse(Lklhd_space_proj, nextLklhdSpace);
%         if isnan(sim)
%             disp('blarg')
%         end
% 
%         
%         corCoef =corrcoef(Lklhd_space_proj,nextLklhdSpace);
%         sim = corCoef(1,2);
%         %[~, sim] = psnr(aa, nextLklhdSpace);
%         sim = max(AScompare(:));   
%         
%         difx = sum(abs(x1dist-x2dist))/length(x1dist);
%         dify = sum(abs(y1dist-y2dist))/length(y1dist);
%         
%         sim = sum(sum(abs(aa-nextLklhdSpace)))/prod(size(nextLklhdSpace));
%         
        simValue(jj)=  sim;
        
%         locs = cat(1,simStruct.hydrophone_struct.location)
%         
%         % For plotting/Debugging/Making nice figures
%         subplot(4,2,(jj*2)); imagesc(simStruct.array_struct.longrid,...
%             simStruct.array_struct.latgrid, nextLklhdSpace); axis xy
%         hold on
%         scatter(locs(:,2), locs(:,1), 'kd', 'filled')
%         scatter(locs(5,2), locs(5,1), 'rd', 'filled')
%         xlabel('Longitude')
%         ylabel('Latitude')
%         
%         subplot(4,2,(jj*2)-1); imagesc(simStruct.array_struct.longrid,...
%             simStruct.array_struct.latgrid, Lklhd_space_proj); axis xy
%         
%         hold on; scatter(locs(:,2), locs(:,1), 'kd', 'filled')
%         scatter(locs(5,2), locs(5,1), 'rd', 'filled')
%         xlabel('Longitude')
%         ylabel('Latitude')
%         title(['Elapsed Time: ', num2str(round(projAmbSurfCurrent.deltaSec(projIdx))), ' sec'])
%         


        
    end

    
    Sim_mat(ii, ii:ii+length(simValue)-1) = simValue;
    Sim_mat(ii:ii+length(simValue)-1,ii) = simValue;
    Sim_mat(ii,ii) = 1;
    
    disp([num2str(ii), ' of ', num2str(length(arrivalArray))])
end







  








end