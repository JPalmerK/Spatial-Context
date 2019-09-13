
function gplDet2RavenSelTable(hyd, fname)




RavenTable = table();
% basic frequency information
f_low =  hyd(5).detection.parm.freq_lo; % frequency limits
f_high = hyd(5).detection.parm.freq_hi;



if ~isfield(hyd(5).detection.calls, 'f_high')
    
    fr=linspace(hyd(5).detection.parm.freq_lo,...
    hyd(5).detection.parm.freq_hi, ...
    hyd(5).detection.parm.bin_hi-...
    hyd(5).detection.parm.bin_lo+1);
    
    
    
    for ii=1:length(hyd)
        
        if ~isempty(hyd(ii).detection)
        for jj=1:length(hyd(ii).detection.calls)
            
            
            cm_max =  GPL_full('cm_max',jj, hyd(ii).detection.calls);
            non_zero_sums = find(sum(cm_max'));
            
            % If there were data present then save the detections
            if~isempty(non_zero_sums)
                
                
                fmin = fr(min(non_zero_sums));
                fmax = fr(max(non_zero_sums));
                
                
                
            else
                fmin =f_low;
                fmax =f_high;
                
                
            end
            
            hyd(ii).detection.calls(jj).flow = fmin;
            hyd(ii).detection.calls(jj).high = fmax;
        end
        end
    end
end


% Make the Raven Table
hyd_idx =1;

for ii =1:length(hyd)
    
if ~isempty(hyd(hyd_idx).detection.calls)
    calls = struct2table(hyd(hyd_idx).detection.calls);
    
    n_calls = height(calls);
    call_ids =height(RavenTable)+1:height(RavenTable)+n_calls;
   
   
    % Start Time
    start_times = calls.start_time/2000;
   
    % End Times (DOOoooOOOoOOOoOOM!)
    end_times = calls.end_time/2000;
   
    % Get the high and low frewuency based on the spectrogram
    % parameters
    low_f =calls.flow;
    high_f = calls.high;
        
   
    Selection =  [height(RavenTable)+1:height(RavenTable)+n_calls]';
    View =repmat({'Spectrogram'},[n_calls,1]);
    Channel =  repmat(ii, [n_calls,1]);
    ClusterId = ones([n_calls, 1]);
    matlabDate = calls.julian_start_time;
    thresh = calls.calltype_1_score;
   
    aa = table(Selection,...
        View, ...
        Channel,...
        start_times,...
        end_times,...
        low_f,...
        high_f,...
        ClusterId,...
        matlabDate,...
        thresh,...
        'VariableNames',{'Selection', 'View', 'Channel',...
        'BeginS', 'EndS', 'LowF', 'HighF', 'ClusterId',...
        'MtlbDtStr', 'Trheshold'});
    

    RavenTable =[RavenTable; aa];
end

hyd_idx =hyd_idx+1;
   
end




writetable(RavenTable, fname, 'Delimiter', '\t', 'WriteVariableNames', false)    

header=strcat('Selection\tView\tChannel\tBegin Time (s)\tEnd Time (s)',...
    '\tLow Freq (Hz)\tHigh Freq (Hz)\',...
    'tClusterID\tMT1\n');


% export the file
% - Read input file data content (without header).
 fid_in = fopen(fname, 'r') ;             % Open input file for reading.
 fgetl(fid_in) ;                                 % Skip header line in input file.
 content = fread(fid_in) ;                       % Read rest of input file.
 fclose(fid_in) ;                                % Close input file.
 
 
 % - Write output file: new header + previous data content.
 fid_out = fopen(fname, 'w') ;  % Open input file for writing.
 fprintf(fid_out, header') ;       % Output new header.
 fwrite(fid_out, content) ;                      % Output previous data content.
 fclose(fid_out) ;  
end
 
 
 
 
 
 
 
 
 
 
 