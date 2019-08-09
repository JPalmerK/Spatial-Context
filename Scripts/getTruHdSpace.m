        function averageLklhd_space = getTruHdSpace(obj, idx, sig_tot)
            
            % Get the observed TDOA space and normalize
            averageLklhd_space = getAvLkHdSpace(obj, idx);
            
            % set up the vectorization for the PDF
            sigma = ones(size(averageLklhd_space)).*sig_tot*sqrt(2*pi);
             
%             % Create ambiguity surface and normalize
%             averageLklhd_space = normpdf(averageLklhd_space, 0, sigma)./...
%                 normpdf(0, 0, sigma);
            
            % Eva comment
            averageLklhd_space = normpdf(averageLklhd_space, 0, sigma).*sigma;
            
            
            if ndims(averageLklhd_space)>1
                
                % sum along third axis, will be normalized later
                %averageLklhd_space = nanmean(averageLklhd_space,3);
                
                temp_ave = averageLklhd_space;
                temp_ave(isnan(temp_ave))=.001;
                % Eva Edit
                averageLklhd_space = prod(temp_ave,3);
                
            end
            
            
            
        end