        function averageLklhd_space = getTruHdSpaceProd(obj, idx, sig_tot)
        % Returns a single grid containing ambiguity surface for a given call deteced
        % by N hydrophones. 
        
            % Get the observed TDOA space and normalize
            averageLklhd_space = getAvLkHdSpace(obj, idx);
            
            % set up the vectorization for the PDF
            sigma = ones(size(averageLklhd_space)).*sig_tot*sqrt(2*pi);
             
%             % Create ambiguity surface and normalize
%             averageLklhd_space = normpdf(averageLklhd_space, 0, sigma)./...
%                 normpdf(0, 0, sigma);
            
            % Eva comment
            averageLklhd_space = normpdf(averageLklhd_space, 0, sigma).*sigma*sqrt(2*pi);
            
            
            if ndims(averageLklhd_space)>1
                
                % sum along third axis, will be normalized later
                %averageLklhd_space = nanmean(averageLklhd_space,3);
                
                temp_ave = averageLklhd_space;
                temp_ave(isnan(temp_ave))=.001;
                
                averageLklhd_space = prod(temp_ave,3, 'omitnan');
                
            end
            
            
            
        end