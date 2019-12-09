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

callsDurFFTStart = zeros(1,length(hyd(local_array(1)).detection.calls));
callsStartEnd = zeros(length(hyd(local_array(1)).detection.calls),2);

for k=1:length(hyd(local_array(1)).detection.calls)
    
    callsDurFFTStart(k) = (hyd(local_array(1)).detection.calls(k).end_time ...
        - hyd(local_array(1)).detection.calls(k).start_time +1)/skip;
    callsStartEnd(k,:)=[hyd(local_array(1)).detection.calls(k).start_time,...
        hyd(local_array(1)).detection.calls(k).end_time];
end

for i1=2:length(local_array)
    clear sf1 sz1
    disp(['Checking hydrophone ', num2str(i1), 'of ', num2str(length(local_array))])
    %i1
    % backwards compatible
    if isfield(array_struct(array_number),'phase') ==0
        %phase=localize_struct.phase;
        phase=localize_struct.parm.phase; % KJP change 10/22/19
    else
        phase=array_struct(array_number).phase(i1-1);
    end
    
    callsSlaveDurfftSamp = zeros(1,length(hyd(local_array(i1)).detection.calls));
    callsSlvStartEnd = zeros(length(hyd(local_array(i1)).detection.calls),2);
    
    for k=1:length(hyd(local_array(i1)).detection.calls)
        
        % duration of call in bins
        callsSlaveDurfftSamp(k) = (hyd(local_array(i1)).detection.calls(k).end_time ...
            - hyd(local_array(i1)).detection.calls(k).start_time +1)/skip;
        
        % start and end call times in samples
        callsSlvStartEnd(k,:)=[hyd(local_array(i1)).detection.calls(k).start_time,...
            hyd(local_array(i1)).detection.calls(k).end_time];
    end
    
    % array length bins (samples/fft skip) for last call
    diss=zeros(1,ceil(callsSlvStartEnd(end,2)/skip));
    
    % Bin index of start times
    Bindex=round((callsSlvStartEnd(:,1)-1)/skip);
    
    % For each call, 
    for k=1:length(hyd(local_array(i1)).detection.calls)
        % energy sum in each time bin
        den=sum(GPL_full('cm_max',k, ...
            hyd(local_array(i1)).detection.calls).^(2*pow));
        
        % Fill in the energy sums?
        diss(Bindex(k)+1:Bindex(k)+callsSlaveDurfftSamp(k))=den;
    end
    % Find calls on the master with start times greater than the first
    % Slave index
    kstart=find(callsStartEnd(:,1)-phase*skip>callsSlvStartEnd(1,1));
    % Set the start index to the first call on the master that occures
    % after the first call on the slave hydrophone
    if(~isempty(kstart > 0))
        kstart=kstart(1);
    end
    % last call on the master worth looking at
    kend=find(callsStartEnd(:,2)+phase*skip<=callsSlvStartEnd(end,2));
    % Index of the last call to consider
    if(~isempty(kend > 0))
        kend=kend(end)-num_calls;
    end
    
    % Array,two times the phase (rows), columns as many as detections
    pprod=zeros(2*phase+1,...
        length(hyd(local_array(1)).detection.calls));
    
    % If there is something in the index starts and end
    if (~isempty(kstart > 0) && ~isempty(kend > 0) && kend >= kstart)
        
        clear num_bs den_b den_s % reset for each master-slave
        
        if num_calls>1
            
            kk=kstart;
            for ic=1:num_calls
                base=GPL_full('cm_max',kk+ic,hyd(local_array(1)).detection.calls).^pow;
                tim=callsStartEnd(kk+ic,:);
                %[b1,b2]=size(base);
                baseSize = size(base,2);
                
                [x,y]=find(base);
                lower=min(x); upper=max(x);
                base=base(lower:upper,min(y):max(y));
                den_base(ic)=sum(sum(base.^2));
                
                interval=[tim(1)-phase*skip,tim(2)+phase*skip];
                k1=find(callsSlvStartEnd(:,1)<=interval(1));k1=k1(end);
                k2=find(callsSlvStartEnd(:,2)>=interval(2));k2=k2(1);
                
                subset=zeros(nfreq,round((callsSlvStartEnd(k2,2)-callsSlvStartEnd(k1,1)+1)/skip) );
                for j=k1:k2
                    start=round((callsSlvStartEnd(j,1)-callsSlvStartEnd(k1,1))/skip); % round in case of noninteger
                    add=GPL_full('cm_max',j,hyd(local_array(i1)).detection.calls).^pow;
                    
                    %[s1,s2]=size(add);
                    s2 = size(add,2);
                    subset(:,start+1:start+s2)=add;
                end
                
                int_close=[callsSlvStartEnd(k1,1),callsSlvStartEnd(k2,2)];
                cut1=round((interval(1)-int_close(1))/skip);
                
                subset=subset(:,cut1+1:cut1+baseSize+2*phase);
                for j=1:(2*phase)+1
                    slate=subset(:,j:baseSize+j-1);
                    slate=slate(lower:upper,min(y):max(y));
                    num_bs(ic,j)=sum(sum(base.*slate));
                end
                
            end % ic loop
            
            den_b=sum(den_base);
            
            den_s = zeros(1,2*phase+1);
            for j=1:(2*phase+1)
                lm=round([(callsStartEnd(kk+1,1)-1)/skip,callsStartEnd(kk+num_calls,2)/skip])+j-1-phase;
                den_s(j)=sum(diss(lm(1):lm(2)));
            end
            pprod(:,kk)=(sum(num_bs)./den_s.^.5/den_b^.5)';
            k=find(den_s==0);
            pprod(k,kk)=0;  % safety valve!
            for kk=kstart+1:kend
                disp([num2str(kk) ' of ' num2str(kend)])
                
                %             if mod(kk,500)==0
                %                 %kk;
                %             end
                
                base=GPL_full('cm_max',kk+num_calls,hyd(local_array(1)).detection.calls).^pow;
                tim=callsStartEnd(kk+num_calls,:);
                %[b1,b2]=size(base);
                baseSize = size(base,2);
                
                [x,y]=find(base);
                lower=min(x); upper=max(x);
                base=base(lower:upper,min(y):max(y));
                % add on one more unit, knock off old
                for ic=1:num_calls-1
                    den_base(ic)=den_base(ic+1);
                end
                den_base(num_calls)=sum(sum(base.^2));
                
                interval=[tim(1)-phase*skip,tim(2)+phase*skip];
                k1=find(callsSlvStartEnd(:,1)<=interval(1));
                if isempty(k1)
                    test = 'Hi';
                end
                k1=k1(end);
                k2=find(callsSlvStartEnd(:,2)>=interval(2));
                if isempty(k2)
                    test = 'Hi';
                end
                k2=k2(1);
                
                subset=zeros(nfreq,round((callsSlvStartEnd(k2,2)...
                    -callsSlvStartEnd(k1,1)+1)/skip)); % round in case of noninteger
                for j=k1:k2
                    start=round((callsSlvStartEnd(j,1)...
                        -callsSlvStartEnd(k1,1))/skip); % round in case of noninteger
                    add=GPL_full('cm_max',j,...
                        hyd(local_array(i1)).detection.calls).^pow;
                    %[s1,s2]=size(add);
                    s2 = size(add,2);
                    subset(:,start+1:start+s2)=add;
                end
                
                
                int_close=[callsSlvStartEnd(k1,1),callsSlvStartEnd(k2,2)];
                cut1=round((interval(1)-int_close(1))/skip);
                subset=subset(lower:upper,cut1+1:cut1+baseSize+2*phase);
                
                num_bs(1:num_calls-1,:)=num_bs(2:num_calls,:);
                
                [b1,baseSize]=size(base);base=reshape(base,1,b1*baseSize);
                for j=1:(2*phase)+1
                    slate=reshape(subset(:,j+min(y)-1:j+max(y)-1),b1*baseSize,1);
                    num_bs(num_calls,j)=base*slate;
                end
                
                % put it all together
                den_b=sum(den_base);
                
                for j=1:(2*phase+1)
                    lm=round([(callsStartEnd(kk+1,1)-1)/skip,callsStartEnd(kk+num_calls,2)/skip])+j-1-phase;
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
                
                % Cleaned spectrogram
                base=GPL_full('cm_max',...
                    kk+1,...
                    hyd(local_array(1)).detection.calls).^pow;
                % start and end time of the call
                tim=callsStartEnd(kk+1,:);
                %[b1,b2]=size(base);
                baseSize = size(base,2);
                
                % Locations where there is energy in the cleaned
                % spectrogram
                [x,y]=find(base);
                
                % Trim the cleaned spectrogram times
                lower=min(x); 
                upper=max(x);
                base=base(lower:upper,min(y):max(y));
                den_base=sum(sum(base.^2));
                
                % Interval for the call start and end
                interval=[tim(1)-phase*skip, tim(2)+phase*skip];
                
                % Find calls on the slave that fall within the start and
                % end times
                k1=find(callsSlvStartEnd(:,1)<=interval(1));
                k1=k1(end);
                k2=find(callsSlvStartEnd(:,2)>=interval(2));
                k2=k2(1);
                
                subset=zeros(nfreq,round((callsSlvStartEnd(k2,2)-callsSlvStartEnd(k1,1)+1)/skip) );
                for j=k1:k2
                    start=round((callsSlvStartEnd(j,1)-callsSlvStartEnd(k1,1))/skip); % round in case of noninteger
                    add=GPL_full('cm_max',j,hyd(local_array(i1)).detection.calls).^pow;
                    
                    %[s1,s2]=size(add);
                    s2 = size(add,2);
                    subset(:,start+1:start+s2)=add;
                end % j loop
                
                int_close=[callsSlvStartEnd(k1,1),callsSlvStartEnd(k2,2)];
                cut1=round((interval(1)-int_close(1))/skip);
                subset=subset(:,cut1+1:cut1+baseSize+2*phase);
                for j=1:(2*phase)+1
                    slate=subset(:,j:baseSize+j-1);
                    slate=slate(lower:upper,min(y):max(y));
                    num_bs(j)=sum(sum(base.*slate));
                end % j loop
                
                for j=1:(2*phase+1)
                    lm=round([(callsStartEnd(kk+1,1)-1)/skip,callsStartEnd(kk+num_calls,2)/skip])+j-1-phase;
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


