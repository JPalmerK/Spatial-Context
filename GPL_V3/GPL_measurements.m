function [GPL_struct] = GPL_measurements(GPL_struct,sp,sp_whiten,quiet_fft, start,finish,parm);
%
%  quiet_fft=2*quiet_fft.^2; %multiplication by 2 corrects for negative frequencies
%  quiet_fft=quiet_fft/parm.sample_freq/sum(hamming(parm.fftl).^2);%normalize by freq bin width (Regina added)
%  quiet_fft=sum(quiet_fft,1);%integrate across frequency band of interest
%  quiet_fft=quiet_fft*parm.sample_freq/parm.fftl;%scale by df (Hz)
%  quiet_fft=sqrt(mean(quiet_fft));%(RMS in uPa) to get dB, take 20*log10(RMS)

for k=1:length(start)
    
    % access & reconstruct cm/cm_max matrices as needed
    
    if(isfield(GPL_struct,'cm')); cm = GPL_full('cm',k,GPL_struct); end
    if(isfield(GPL_struct,'cm_max')); cm_max =  GPL_full('cm_max',k,GPL_struct); end
    if(isfield(GPL_struct,'cm_max2')); cm_max2 = GPL_full('cm_max2',k,GPL_struct); end
    
    
    
    
    if(parm.measure.slope == 1 || parm.measure.cm_max_duration ==1)
        nonz=find(sum(cm_max));
        kk=length(nonz);
        
        GPL_struct(k).cm_max_duration_time = kk * (parm.skip/parm.sample_freq);
        GPL_struct(k).cm_max_duration_bin = kk ;
    end
    
    if(parm.measure.slope == 1 || parm.measure.cm_max2_duration ==1)
        nonz=find(sum(cm_max2));
        kk=length(nonz);
        GPL_struct(k).cm_max2_duration_time = kk * (parm.skip/parm.sample_freq);
        GPL_struct(k).cm_max2_duration_bin = kk ;
    end
    
    if(parm.measure.slope==1)
        GPL_struct(k).slope_ci=nan;
        GPL_struct(k).slope=nan;
        
        nonz=find(sum(cm_max));
        kk=length(nonz);
        
        if kk>2
            
            freq=[0:parm.fftl/2]/parm.fftl*parm.sample_freq;
            freq=freq(parm.bin_lo:parm.bin_hi)';
            
            wt=sum(cm_max');  kf=find(wt);
            kf=kf([1,end]);freq_range=freq(kf);
            
            wt=sum(cm_max');fnd=(wt*freq)/sum(wt);
            [ys,xs,ws]=find(cm_max);
            
            [fitresult, gof] = weighted_slope(xs, ys, ws.^2);
            ci = confint(fitresult,0.95);
            confidence = abs(ci(1)-ci(2));
            
            GPL_struct(k).slope_ci=confidence;
            GPL_struct(k).slope = fitresult.p1;
            
            
        end
    end
    
    if(parm.measure.cm_max_bandwidth == 1)
        %%% Record the number of frequency bins in contour
        jj=length(find(sum(cm_max')));
        GPL_struct(k).cm_max_bandwidth_freq = jj*(parm.sample_freq/parm.fftl);
        GPL_struct(k).cm_max_bandwidth_bin = jj;
    end
    
    if(parm.measure.cm_max2_bandwidth == 1)
        jj=length(find(sum(cm_max2')));
        GPL_struct(k).cm_max2_bandwidth_freq = jj*(parm.sample_freq/parm.fftl);
        GPL_struct(k).cm_max2_bandwidth_bin = jj;
    end
    
    
    
    if(parm.measure.spec_noise==1)
        dt=parm.skip/parm.sample_freq;%dt=(nfft-noverlap)/fs
        timebins=ceil(5/dt);%how many time bins are in 5 seconds?
        skipbins=ceil(.5/dt);
        if start(k)-skipbins-timebins>=1
            sp_noise=sp(:,((start(k)-skipbins)-timebins):(start(k)-skipbins));%skip a bin between call and noise, take 2 seconds of time
        elseif start(k)-skipbins>=1
            sp_noise=sp(:,1:(start(k)-skipbins));%skip bins between call and noise, take 2 seconds of time
            [~,colsp]=size(sp_noise);
            sp_noise=[sp_noise,sp(:,finish(k)+skipbins:finish(k)+skipbins+timebins-colsp)];%add some noise from after the call
        else
            sp_noise=sp(:,finish(k)+skipbins:finish(k)+skipbins+timebins);%use time after call if call is at start of segment
        end
            sp_noise=2*sp_noise.^2; %multiplication by 2 corrects for negative frequencies
        sp_noise=sp_noise/parm.sample_freq/sum(hamming(parm.fftl).^2);%normalize by freq bin width (Regina added)
        sp_noise=sum(sp_noise,1);%integrate across frequency band of interest
        sp_noise=sp_noise*parm.sample_freq/parm.fftl;%scale by df (Hz)
        sp_noise=sqrt(mean(sp_noise));%(RMS in uPa) to get dB, take 20*log10(RMS)
                
        GPL_struct(k).spec_noise = sp_noise;
    end
    
    if(parm.measure.spec_rl==1)
        
        if length(find(cm)) > 5
            
            sp_subset=sp(:,start(k):finish(k));
            
            [spec_rl]=estimate_rl(cm,sp_subset,parm);
        else
            spec_rl=nan;
        end
                
        GPL_struct(k).spec_rl = spec_rl;
    end
    
    
end

%GPL_struct=rmfield(GPL_struct,'cm');



