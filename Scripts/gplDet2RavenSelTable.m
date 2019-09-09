
function gplDet2RavenSelTable(hyd, fname)




RavenTable = table();

% frequency information
f_low =  hyd(1).detection.parm.freq_lo; % frequency limits
f_high = hyd(1).detection.parm.freq_hi;
df = (f_high-f_low)/(hyd(1).detection.parm.bin_hi-hyd(1).detection.parm.bin_lo);


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
   
    aa = table(Selection,...
        View, ...
        Channel,...
        start_times,...
        end_times,...
        low_f,...
        high_f,...
        ClusterId,...
        matlabDate,...
        'VariableNames',{'Selection', 'View', 'Channel',...
        'BeginS', 'EndS', 'LowF', 'HighF', 'ClusterId','MtlbDtStr'});
    

    RavenTable =[RavenTable; aa];
end

hyd_idx =hyd_idx+1;
   
end




writetable(RavenTable, fname, 'Delimiter', '\t', 'WriteVariableNames', false)    

header=strcat('Selection\tView\tChannel\tBegin Time (s)\tEnd Time (s)',...
    '\tLow Freq (Hz)\tHigh Freq (Hz)\tClusterId\tMtLbDatestr\n');


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
 
 
 
 
 
 
 
 
 
 
 