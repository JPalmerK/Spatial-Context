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
    
%     disp(['End section ',...
%         num2str(round(((j-1)*parm.nrec+parm.nrec)/2000/60,2)),...
%         ' min'])
%     
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
           
            
                fmin = fr(min(non_zero_sums));
                fmax = fr(max(non_zero_sums));

        
            
        else 
            fmin = parm.sum_freq_lo;
            fmax = parm.sum_freq_hi;
        
       
        end
       
        GPL_struct(k).flow = fmin;
        
        GPL_struct(k).high = fmax;
        
    end
    
      calls=[calls,GPL_struct];
    
    
end


end
