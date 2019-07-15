% Create clusters of data from GPL localization struct
classdef clusterClassGPL <simulationClass
    properties

         
        limitTime % to limit the time over which to analyse the data fill here (seconds)
        
        % Exporting and comparing
        ravenColumnHeadings ={'Selection', 'View', 'Channel', 'Begin Time (s)', 'End Time (s)', 'Low Freq (Hz)', 'High Freq (Hz)'};
        RavenTable
        
        % Importing data (GPL or selection table)
        localize_struct % Localize struct if using GPL
        raventLoc % Location of Raven selection table
        
        % For referencing time of GPL detections
        calls
        callParm
        
    end
    
    
    methods
        function UpdateArrTable(obj)
            % Create the arrival table based on either GPL output or raven
            % selection table
            
            if isempty(obj.localize_struct) & isempty(obj.raventLoc)
                disp('Raven Selection Table or GPL localization Structure needed')
                return
            
            elseif  isempty(obj.localize_struct) & ~isempty(obj.raventLoc)
                updateArrTableRaven(obj) % Not implemented
                disp('Ravent table access not yet implemented')
            
            elseif ~isempty(obj.localize_struct) & isempty(obj.raventLoc)
                updateArrTableGPL(obj)
            
            else ~isempty(obj.localize_struct) & ~isempty(obj.raventLoc)
                disp('Raven selection table OR GPL localize struct needed')
                return
            end
        end
        
        %%  Create the array table from Ravent Textfile
        function updateArrTableRaven(obj)
            
            ravenTable = readtable(obj.ravenLoc)
        
        
        end
        
        
        %% Create a Raven Compatable Selection table 
        function RavenTable = exportRavenTxt(obj, outputDir)
            % Function ofr creating and exporting detections as Raven
            % compatable .txt file
            

            parent = obj.array_struct.master;
            arivalArr = obj.arrivalArray;
            
            RavenTable = array2table(zeros(0,8),...
                'VariableNames',{'Selection', 'View', 'Channel',...
                'BeginS', 'EndS', 'LowF', 'HighF','ClusterId'});
            
            % Loop through the child indexes and create the arrival tables
            hyds = [obj.array_struct.master, obj.array_struct.slave(obj.child_idx)];
            
            % Calls 
            calls = struct2table(obj.calls);
            
            % frequency information
            f_low = obj.callParm.freq_lo; % frequency limits
            f_high = obj.callParm.freq_hi;
            
            df = (f_high-f_low)/(obj.callParm.bin_hi-obj.callParm.bin_lo);
            
           
            
            for ii =1:length(hyds)
                n_calls =sum(~isnan(arivalArr(:,ii)));
                call_ids =obj.localize_struct.hyd(parent).dex(logical(~isnan(arivalArr(:,ii))));
                
                % Error in code
                call_ids = call_ids;
                
                % Start Time
                start_times = calls.start_time(call_ids)/obj.fs;
                
                % End Times (DOOoooOOOoOOOoOOM!)
                end_times = calls.end_time(call_ids)/obj.fs;
                
                % Get the high and low frewuency based on the spectrogram
                % parameters
                low_f = zeros(length(end_times),1);
                high_f = low_f;
                
                for jj=1:n_calls
                    [rr, cc]=ind2sub(calls.cm(call_ids(jj)).size, calls.cm(call_ids(jj)).index);
                    low_f(jj) = f_low+(df*min(rr));
                    high_f(jj) = f_low+(df*max(rr));
                end
                
                
                Selection =  [height(RavenTable)+1:height(RavenTable)+n_calls]';
                View =repmat({'Spectrogram'},[n_calls,1]);
                Channel =  repmat(hyds(ii), [n_calls,1]);
                ClusterId = obj.Cluster_id(logical(~isnan(arivalArr(:,ii))));
                
                aa = table(Selection,...
                    View, ...
                    Channel,...
                    start_times,...
                    end_times,...
                    low_f,...
                    high_f,...
                    ClusterId,...
                    'VariableNames',{'Selection', 'View', 'Channel',...
                    'BeginS', 'EndS', 'LowF', 'HighF', 'ClusterId'});
                
                RavenTable =[RavenTable; aa];
                
            end
            
            
            obj.RavenTable =RavenTable;
            
            if nargin==2
                fnameloc = outputDir;
            end
            
        end
        
        %% Create the array table from GPL data
        function updateArrTableGPL(obj)
            % Update the arrival table using GPL detections
            
            % Parent hydrophone
            parent = obj.array_struct.master;
            child_hyd = obj.array_struct.slave;
            
            ParentArrival = obj.localize_struct.hyd(parent).rtimes'./obj.fs;
            TDOA =obj.localize_struct.hyd(parent).delays;
            
            ArrivalSec =zeros(length(TDOA), (length(obj.hydrophone_struct)))/0;
            ArrivalSec(:,parent) =ParentArrival;
            CrossScores =zeros(length(TDOA), (length(obj.hydrophone_struct)))/0;
            
            for ii=1:length(child_hyd)
                ArrivalSec(:,child_hyd(ii)) = (ParentArrival+TDOA(:,ii)).*...
                    (obj.localize_struct.hyd(parent).cross_score(:,ii)./...
                    obj.localize_struct.hyd(parent).cross_score(:,ii));
                CrossScores(:,child_hyd(ii)) = obj.localize_struct.hyd(parent).cross_score(:,ii);
            end
            
            
            x=squeeze(obj.localize_struct.hyd(parent).coordinates(2,1,:));
            x(:,2)=squeeze(obj.localize_struct.hyd(parent).coordinates(2,2,:));
            
            % Index of the detection
            idx = obj.localize_struct.hyd(parent).dex';
            
            %
            % Create the arrivals table (at) from the localize structure
            at = table(...
                ArrivalSec,...
                x,...
                CrossScores,...
                TDOA,...
                idx,...
                'VariableNames',{'ArrivalSec', 'Location','CrossScore', 'TDOA', 'dex'});
            at.ID = zeros(height(at),1)/0;
            % Clear out any detections greater than Xkm from the receiver

            
            
            % limit the amount of time considered
            if ~isempty(obj.limitTime)
                good_idx = logical(cumsum(diff(ArrivalSec(:, parent)))<obj.limitTime);
                at = at(good_idx,:);
            end
            
            obj.arrivalTable =at;
        end
        
        %% Draw the locations (where knonw) and color by cluster ID
        function drawAgents(obj)
            
            % Update the arrival array
            if isempty(obj.arrivalArray)
                UpdateArrArray(obj);
            end
            
            temp_array = obj.arrivalArray;
            temp_array(:,end) =obj.Cluster_id;
            % Get rid of rows without a location
            good_idx = logical(~isnan(temp_array(:,end-1)));
            temp_array = temp_array(good_idx,:);
            
            
            duration = temp_array(end,1)-temp_array(1,1);
            newTimescale = round(temp_array(:,1)-temp_array(1,1))+1;
            
            
            
            % Hydrophone locations
            hyd_table = struct2table(obj.hydrophone_struct);
            TimeColorVals = parula(duration+5);
            ColorVals = lines(max(temp_array(:,end)));
            
            child_hyds = obj.array_struct.slave(obj.child_idx)
            
            
            figure;
            subplot(2,1,1)
            hold on
            scatter(temp_array(:,end-1),...
                temp_array(:,end-2),[],...
                TimeColorVals(newTimescale,:), 'f')
            scatter(hyd_table.location(:,2),...
                hyd_table.location(:,1), 80,...
                'k', 'filled', 'd')
            scatter(...
                hyd_table.location([obj.array_struct.master,...
                obj.array_struct.slave(obj.child_idx)],2),...
                hyd_table.location([obj.array_struct.master,...
                obj.array_struct.slave(obj.child_idx)],1),...
                'r', 'filled', 'd')
            
            title('TOA on Parent')
            
            subplot(2,1,2)
            hold on
            scatter(temp_array(:,end-1),...
                temp_array(:,end-2),[],...
                ColorVals(temp_array(:,end),:), 'f')
            scatter(hyd_table.location(:,2),...
                hyd_table.location(:,1), 80,...
                'k', 'filled', 'd')
            scatter(...
                hyd_table.location([obj.array_struct.master,...
                obj.array_struct.slave(obj.child_idx)],2),...
                hyd_table.location([obj.array_struct.master,...
                obj.array_struct.slave(obj.child_idx)],1),...
                'r', 'filled', 'd')
            
            
            
            
            
            
            
        end
    end
end