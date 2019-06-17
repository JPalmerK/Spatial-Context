function [calls]=process_RW_sub(data, parm,start_date, num_slates_file); 
     
    tic 
    P=parm.sample_freq; 
     
     
    scheck=length(data); 
    if(scheck < P*1*10) 
        data=[data,zeros(1,(P*1*10)-scheck)]; 
    end 
     
    data=data'; 
     
    calls=[]; 
    for(j=1:num_slates_file)  %change this depending on length of record reading in per file 
        sub_data=data(1+(j-1)*parm.nrec:((j-1)*parm.nrec)+parm.nrec); 
         
         
        [GPL_struct]=GPL_v3(sub_data,parm); 
         
         
        for(k=1:length(GPL_struct)) 
             
            GPL_struct(k).start_time=GPL_struct(k).start_time+((j-1)*parm.nrec); 
            GPL_struct(k).end_time=GPL_struct(k).end_time+((j-1)*parm.nrec); 
             
            GPL_struct(k).julian_start_time=start_date+datenum(0,0,0,0,0,GPL_struct(k).start_time/parm.sample_freq); 
            GPL_struct(k).julian_end_time=start_date+datenum(0,0,0,0,0,GPL_struct(k).end_time/parm.sample_freq); 
             
        end 
         
        calls=[calls,GPL_struct]; 
         
    end 
     