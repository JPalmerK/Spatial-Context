load call1; load call2; load call3; load call4; load call5; load call6; load call7; load call8; load call9; load call10; load call11;
%load case2;


a(1).call=call1;
a(2).call=call2;
a(3).call=call3;
a(4).call=call4;
a(5).call=call5;
a(6).call=call6;
a(7).call=call7;
a(8).call=call8;
a(9).call=call9;
a(10).call=call10;

cc=hsv(10);

for(j=1:100)
    
    begin = round(random('uniform',1,length(case1)-(60*10000)))
    noise=case1(begin+1:begin+60*10000);
    noise=noise-mean(noise);
    for(jj=1:10)
        
        ucall=a(jj).call;
        ucall=ucall./max(ucall);
        call=ucall;
        
        
        
        
        
        call=[zeros(1,20*10000),call',zeros(1,40*10000-length(call))];
        call=call';
        call=call./(max(call));
        db_extract=[];
        db_orig=[];
        db_extract_pp=[];
        db_orig_pp=[];
        db=[];
        
        
        for(factor=[5,4,3,2,1,.9,.8,.7,.6,.5,.4,.3,.2,.1,.05,.01, .001, .0001])
            factor
            combined=noise+(call*(factor)*max(noise));
            
            %sp=spmaker(combined,10000,1000,2048,1);
            
            
            [times,dspec,sndtrk,baseline0,diagnostic,SNR,clean_call]=newGPL_iter2_smee(combined,parm);
            
            
            
            if(isempty(clean_call)==0)
                %figure(21); imagesc(clean_call(1).spec)
                
                temp=[];
                for(kk=1:length(clean_call))
                    temp=[temp;clean_call(kk).times];
                end
                
                [k1 k2]=min(abs(temp(:,1)-200000));
                
                if(k1/10000 < 1)
                    
                    %  figure(1)
                    %  imagesc(20*log10(clean_call(k2).spec));
                    %  drawnow;
                    
                    db_extract=[db_extract;20*log10(clean_call(k2).spec_rl)];
                    
                    test=spmaker(ucall*(factor)*max(noise),10000,150,1000,2048,0,512);
                    test = sqrt(mean(sum(test,1)));
                    db_orig=[db_orig;20*log10(test)];
                    
                    %If doing it from time series
                    %db_orig=[db_orig;20*log10(sqrt(mean(abs(ucall*(factor)*max(noise)).^2)))];
                    
                    
                    db=[db;20*log10(factor)];
                    
                    
                end
            end
        end
        
        
        
        
        %
        % figure; plot(db,100+db_orig,'b');
        % hold on; plot(db,100+db_extract,'g');
        %
        %
        % figure; plot(db,100+db_orig_pp,'b');
        % hold on; plot(db,100+db_extract_pp,'g');
        
        
        figure(22); plot(db,db_orig-db_extract,'color',cc(jj,:));
        drawnow;
        %figure(23); plot(db,db_orig_pp-db_extract_pp,'color',cc(jj,:));
        hold on;
        
    end
end
