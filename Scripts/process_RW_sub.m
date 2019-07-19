function [calls]=process_RW_sub(data, parm,start_date, num_slates_file);

tic
P=parm.sample_freq;


scheck=length(data);
if(scheck < P*1*10)
    data=[data,zeros(1,(P*1*10)-scheck)];
end

data=data';
fr=linspace(parm.freq_lo,parm.freq_hi,parm.bin_hi-parm.bin_lo+1);
calls=[];
for(j=1:num_slates_file)  %change this depending on length of record reading in per file
    sub_data=data(1+(j-1)*parm.nrec:((j-1)*parm.nrec)+parm.nrec);
    
    
    [GPL_struct]=GPL_v3(sub_data,parm);
    
    
    for(k=1:length(GPL_struct))
        
        GPL_struct(k).start_time=GPL_struct(k).start_time+((j-1)*parm.nrec);
        GPL_struct(k).end_time=GPL_struct(k).end_time+((j-1)*parm.nrec);
        
        GPL_struct(k).julian_start_time=start_date+datenum(0,0,0,0,0,GPL_struct(k).start_time/parm.sample_freq);
        GPL_struct(k).julian_end_time=start_date+datenum(0,0,0,0,0,GPL_struct(k).end_time/parm.sample_freq);
        
        
        cm_max =  GPL_full('cm_max',k,GPL_struct);
        non_zero_sums = find(sum(cm_max'));
        
        % If there were data present then save the detections
        if~isempty(non_zero_sums)
        gaps = diff(non_zero_sums);
        first_gap =0;
        
        if sum(gaps)>length(gaps)
            % get the first gap for f min
            first_gap = find(gaps>1,1);
            fmin = fr(first_gap+1);
        else
            fmin = fr(min(non_zero_sums));
        end
        
        
        % Get the second gap for f max
        first_gap = non_zero_sums(first_gap+1:end);
        gaps = diff(first_gap);
        if sum(gaps)> length(gaps)
            last_gap =  max(find(gaps>5));
            fmax = fr(first_gap(last_gap));
        else
            fmax = fr(max(first_gap));
        end
        
        else
            fmin = 0/0;
            fmax = fmin;
        end
        
        GPL_struct(k).flow = fmin;
        
        GPL_struct(k).high = fmax;
        
        
    end
    
    
    
    
end

calls=[calls,GPL_struct];

end
