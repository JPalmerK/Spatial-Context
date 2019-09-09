function [localize_struct] = localize_cross_corr_new(array_struct,hyd,localize_struct,array_number);


% Number of calls to cross correlate
num_calls = localize_struct.parm.num_calls;
pow = localize_struct.parm.pow;

%fftl= hyd(array_struct(array_number).master).detection.parm.fftl;
skip= hyd(array_struct(array_number).master).detection.parm.skip;
nfreq= hyd(array_struct(array_number).master).detection.parm.nfreq;

local_array=[array_struct(array_number).master,...
    array_struct(array_number).slave];

for n=local_array
    if (length(hyd(n).detection.calls)<=num_calls)
        localize_struct.hyd(local_array(1)).cc_matrix = [];
        return;
    end
end

sz2 = zeros(1,length(hyd(local_array(1)).detection.calls));
sf2 = zeros(length(hyd(local_array(1)).detection.calls),2);

for k=1:length(hyd(local_array(1)).detection.calls)
    
    sz2(k) = (hyd(local_array(1)).detection.calls(k).end_time ...
        - hyd(local_array(1)).detection.calls(k).start_time +1)/skip;
    sf2(k,:)=[hyd(local_array(1)).detection.calls(k).start_time,...
        hyd(local_array(1)).detection.calls(k).end_time];
end

for i1=2:length(local_array)
    clear sf1 sz1
    disp(num2str(i1))
    %i1
    % backwards compatible
    if isfield(array_struct(array_number),'phase') ==0
        phase=localize_struct.phase;
    else
        phase=array_struct(array_number).phase(i1-1);
    end
    
    sz1 = zeros(1,length(hyd(local_array(i1)).detection.calls));
    sf1 = zeros(length(hyd(local_array(i1)).detection.calls),2);
    
    for k=1:length(hyd(local_array(i1)).detection.calls)
        
        % duration of call in bins
        sz1(k) = (hyd(local_array(i1)).detection.calls(k).end_time ...
            - hyd(local_array(i1)).detection.calls(k).start_time +1)/skip;
        
        % start and end call times in samples
        sf1(k,:)=[hyd(local_array(i1)).detection.calls(k).start_time,...
            hyd(local_array(i1)).detection.calls(k).end_time];
    end
    
    % array length bins (samples/fft skip) for last call
    diss=zeros(1,ceil(sf1(end,2)/skip));
    
    % Bin index of start times
    dex=round((sf1(:,1)-1)/skip);
    
    % For each call, 
    for k=1:length(hyd(local_array(i1)).detection.calls)
        % energy sum in each time bin
        den=sum(GPL_full('cm_max',k, ...
            hyd(local_array(i1)).detection.calls).^(2*pow));
        
        % Fill in the energy sums?
        diss(dex(k)+1:dex(k)+sz1(k))=den;
    end
    
    kstart=find(sf2(:,1)-phase*skip>sf1(1,1));
    
    if(~isempty(kstart > 0))
        kstart=kstart(1);
    end
    
    kend=find(sf2(:,2)+phase*skip<=sf1(end,2));
    
    if(~isempty(kend > 0))
        kend=kend(end)-num_calls;
    end
    
    
    pprod=zeros(2*phase+1,length(hyd(local_array(1)).detection.calls));
    
    if (~isempty(kstart > 0) && ~isempty(kend > 0) && kend >= kstart)
        
        clear num_bs den_b den_s % reset for each master-slave
        
        if num_calls>1
            
            kk=kstart;
            for ic=1:num_calls
                base=GPL_full('cm_max',kk+ic,hyd(local_array(1)).detection.calls).^pow;
                tim=sf2(kk+ic,:);
                %[b1,b2]=size(base);
                b2 = size(base,2);
                
                [x,y]=find(base);
                lower=min(x); upper=max(x);
                base=base(lower:upper,min(y):max(y));
                den_base(ic)=sum(sum(base.^2));
                
                interval=[tim(1)-phase*skip,tim(2)+phase*skip];
                k1=find(sf1(:,1)<=interval(1));k1=k1(end);
                k2=find(sf1(:,2)>=interval(2));k2=k2(1);
                
                subset=zeros(nfreq,round((sf1(k2,2)-sf1(k1,1)+1)/skip) );
                for j=k1:k2
                    start=round((sf1(j,1)-sf1(k1,1))/skip); % round in case of noninteger
                    add=GPL_full('cm_max',j,hyd(local_array(i1)).detection.calls).^pow;
                    
                    %[s1,s2]=size(add);
                    s2 = size(add,2);
                    subset(:,start+1:start+s2)=add;
                end
                
                int_close=[sf1(k1,1),sf1(k2,2)];
                cut1=round((interval(1)-int_close(1))/skip);
                
                subset=subset(:,cut1+1:cut1+b2+2*phase);
                for j=1:(2*phase)+1
                    slate=subset(:,j:b2+j-1);
                    slate=slate(lower:upper,min(y):max(y));
                    num_bs(ic,j)=sum(sum(base.*slate));
                end
                
            end % ic loop
            
            den_b=sum(den_base);
            
            den_s = zeros(1,2*phase+1);
            for j=1:(2*phase+1)
                lm=round([(sf2(kk+1,1)-1)/skip,sf2(kk+num_calls,2)/skip])+j-1-phase;
                den_s(j)=sum(diss(lm(1):lm(2)));
            end
            pprod(:,kk)=(sum(num_bs)./den_s.^.5/den_b^.5)';
            k=find(den_s==0);pprod(k,kk)=0;  % safety valve!
            for kk=kstart+1:kend
                disp([num2str(kk) ' of ' num2str(kend)])
                
                %             if mod(kk,500)==0
                %                 %kk;
                %             end
                
                base=GPL_full('cm_max',kk+num_calls,hyd(local_array(1)).detection.calls).^pow;
                tim=sf2(kk+num_calls,:);
                %[b1,b2]=size(base);
                b2 = size(base,2);
                
                [x,y]=find(base);
                lower=min(x); upper=max(x);
                base=base(lower:upper,min(y):max(y));
                % add on one more unit, knock off old
                for ic=1:num_calls-1
                    den_base(ic)=den_base(ic+1);
                end
                den_base(num_calls)=sum(sum(base.^2));
                
                interval=[tim(1)-phase*skip,tim(2)+phase*skip];
                k1=find(sf1(:,1)<=interval(1));
                if isempty(k1)
                    test = 'Hi';
                end
                k1=k1(end);
                k2=find(sf1(:,2)>=interval(2));
                if isempty(k2)
                    test = 'Hi';
                end
                k2=k2(1);
                
                subset=zeros(nfreq,round((sf1(k2,2)-sf1(k1,1)+1)/skip)); % round in case of noninteger
                for j=k1:k2
                    start=round((sf1(j,1)-sf1(k1,1))/skip); % round in case of noninteger
                    add=GPL_full('cm_max',j,...
                        hyd(local_array(i1)).detection.calls).^pow;
                    %[s1,s2]=size(add);
                    s2 = size(add,2);
                    subset(:,start+1:start+s2)=add;
                end
                
                
                int_close=[sf1(k1,1),sf1(k2,2)];
                cut1=round((interval(1)-int_close(1))/skip);
                subset=subset(lower:upper,cut1+1:cut1+b2+2*phase);
                
                num_bs(1:num_calls-1,:)=num_bs(2:num_calls,:);
                
                [b1,b2]=size(base);base=reshape(base,1,b1*b2);
                for j=1:(2*phase)+1
                    slate=reshape(subset(:,j+min(y)-1:j+max(y)-1),b1*b2,1);
                    num_bs(num_calls,j)=base*slate;
                end
                
                % put it all together
                den_b=sum(den_base);
                
                for j=1:(2*phase+1)
                    lm=round([(sf2(kk+1,1)-1)/skip,sf2(kk+num_calls,2)/skip])+j-1-phase;
                    den_s(j)=sum(diss(lm(1):lm(2)));
                end
                pprod(:,kk)=(sum(num_bs)./den_s.^.5./den_b^.5)';
                k=find(den_s==0);pprod(k,kk)=0;  % safety valve! solves NaN problem with single calls
                
            end % end of the kk loop
            
        else
            
            %%%%
            for kk=kstart:kend
                
                %             if mod(kk,500)==0
                %                 kk
                %             end
                
                base=GPL_full('cm_max',kk+1,hyd(local_array(1)).detection.calls).^pow;
                tim=sf2(kk+1,:);
                %[b1,b2]=size(base);
                b2 = size(base,2);
                
                [x,y]=find(base);
                lower=min(x); upper=max(x);
                base=base(lower:upper,min(y):max(y));
                den_base=sum(sum(base.^2));
                
                interval=[tim(1)-phase*skip,tim(2)+phase*skip];
                k1=find(sf1(:,1)<=interval(1));k1=k1(end);
                k2=find(sf1(:,2)>=interval(2));k2=k2(1);
                
                subset=zeros(nfreq,round((sf1(k2,2)-sf1(k1,1)+1)/skip) );
                for j=k1:k2
                    start=round((sf1(j,1)-sf1(k1,1))/skip); % round in case of noninteger
                    add=GPL_full('cm_max',j,hyd(local_array(i1)).detection.calls).^pow;
                    
                    %[s1,s2]=size(add);
                    s2 = size(add,2);
                    subset(:,start+1:start+s2)=add;
                end % j loop
                
                int_close=[sf1(k1,1),sf1(k2,2)];
                cut1=round((interval(1)-int_close(1))/skip);
                subset=subset(:,cut1+1:cut1+b2+2*phase);
                for j=1:(2*phase)+1
                    slate=subset(:,j:b2+j-1);
                    slate=slate(lower:upper,min(y):max(y));
                    num_bs(j)=sum(sum(base.*slate));
                end % j loop
                
                for j=1:(2*phase+1)
                    lm=round([(sf2(kk+1,1)-1)/skip,sf2(kk+num_calls,2)/skip])+j-1-phase;
                    den_s(j)=sum(diss(lm(1):lm(2)));
                end % j loop
                
                pprod(:,kk)=(num_bs./den_s.^.5/den_base^.5)';
                k=find(den_s==0);pprod(k,kk)=0;  % safety valve!
                
            end % end of the kk loop
            
        end % which num_calls loop to use
        
    end
    
    %low_index=find(pprod < .1);
    %pprod(low_index)=0;
    pprod(pprod < .1) = 0;
    
    %localize_struct.hyd(local_array(1)).cc_matrix{i1-1}=sparse(pprod);
    %Removed due to indexing error on 1/2/2019
    
    pprod_shift=0*pprod;
    pprod_shift(:,2:end)=pprod(:,1:end-1);
    localize_struct.hyd(local_array(1)).cc_matrix{i1-1}=sparse(pprod_shift);
    
    localize_struct.parm.index_check=1;
    
    %whos pprod
    %local_array(1)
    %i1-1
end % i1 loop


