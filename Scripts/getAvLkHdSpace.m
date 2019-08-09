       function averageLklhd_space = getAvLkHdSpace(obj, callIdx)
            % Returns an averaged likelihood space for a given TDOA
            
            % Pick the set of call delays (number of calls)
            delays = (obj.TDOA_vals(callIdx, :));
            
            
            child_idx = obj.child_idx(~isnan(delays));
            
            averageLklhd_space = zeros([length(obj.array_struct.latgrid),...
                length(obj.array_struct.longrid),...
                length(child_idx)]);
            
            child_hyds = obj.array_struct.slave(obj.child_idx);
            child_hyds = child_hyds(~isnan(delays));
              
            
            for ii=1:length(child_idx)
                
                % Get the expected delay spaces for the hydrophone pairs
                toa_space = (...
                    cell2mat(obj.array_struct.toa_diff(child_idx(ii)+1)));
                
                % Delta TOA space
                averageLklhd_space(:,:,ii) = (toa_space - delays(ii));
                
%                 
%                 figure;
%                 contourf(obj.array_struct.longrid, obj.array_struct.latgrid, (averageLklhd_space(:,:,ii)))
%                 colorbar;
%                 caxis([-3 3])
%                 hold on
%                 scatter(obj.arrivalArray(callIdx,end-1), obj.arrivalArray(callIdx,end-2))
            end

        end

        