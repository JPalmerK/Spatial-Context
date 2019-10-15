function [localize_struct] = localize_LSQ_2D(array_struct, hydrophone_struct,hyd,localize_struct,array_number);
    
    %cutoff_cc=localize_struct.parm.cc_cutoff; % currently unused
    
    local_array=[array_struct(array_number).master, ...
        array_struct(array_number).slave];
    
    if isfield(localize_struct.parm,'min_pairs')
        min_pairs=localize_struct.parm.min_pairs;
    else
        min_pairs=length(local_array)-1;end
    
    % Number of calls to cross correlate
    num_calls=localize_struct.parm.num_calls;
    
    % Number of maxima to cross correlate
    nCallsXcorr=localize_struct.parm.number_maxima;
    
    % LSQ cuttoff
    cutoff=localize_struct.parm.lsq_cutoff;
    
    % Sound speed profile
    ssp=localize_struct.parm.ssp;
    svp=ssp(:,2);
    zz = ssp(:,1);
    
    % Sample parameters
    num_freq=hyd(local_array(1)).detection.parm.nfreq;
    sample_freq=hyd(local_array(1)).detection.parm.sample_freq;
    skip=hyd(local_array(1)).detection.parm.skip;
    
    % Latlon space and E (?!?!)
    E = [6378137.0 0.081819190842613];
    latgrid= array_struct(array_number).latgrid;
    longrid= array_struct(array_number).longrid;
    
    % assign null outputs initially
    score=[]; %
    cross_score=[];
    cc_matrix=[];
    coordinates =[];
    score=[];
    delays=[];
    times=[];
    dex=[];
    cross_score=[];
    rtimes=[];
    
    
    % Get lat/lon and depth of the hydrophone array
    for k=1:length(local_array);
        array(k,:)=hydrophone_struct(local_array(k)).location;
        dep(k)=hydrophone_struct(local_array(k)).depth;
    end
    
    
    
    % Create/pull (??) cross correlation structure for the master vs
    % each slave hydrophone
    for k=1:length(local_array)-1;
        cc_matrix.slv{k}=localize_struct.hyd(local_array(1)).cc_matrix{k};
    end
    
    %(??)
    [sz1,nOtherCalls]=size(cc_matrix.slv{1});
    
    
    % Create empty matrix for delay, rdelay (?) and peak
    delay=zeros(nCallsXcorr,nOtherCalls,length(cc_matrix.slv));
    rdelay=delay;
    peak=delay+nan;
    
    
    % Iterate through the slave hydrophones
    for(jj=1:length(cc_matrix.slv))
        
        % Get the size of the cross correlation matrix
        CC_matSlave=cc_matrix.slv{jj};
        [sz1,nOtherCalls]=size(cc_matrix.slv{jj});
        
        
        
        % Check if phase is present, if not add it
        if isfield(array_struct(array_number),'phase') ==0
            
            % trim if old fashioned
            trm=vdist(hydrophone_struct(local_array(1)).location(1),...
                hydrophone_struct(local_array(1)).location(2),...
                hydrophone_struct(local_array(jj+1)).location(1),...
                hydrophone_struct(local_array(jj+1)).location(2));
            
            trm=ceil(trm/1500*sample_freq/skip);pha(jj)=trm;
            crop=localize_struct.parm.phase-trm;
            
            if crop>0
                CC_matSlave=CC_matSlave(crop+1:end-crop,:);end
            
        else
            pha(jj)=array_struct(array_number).phase(jj);
        end
        
        
        for k=1:nOtherCalls
            ks=find((CC_matSlave(3:end,k)-CC_matSlave(2:end-1,k)).*...
                (CC_matSlave(2:end-1,k)-CC_matSlave(1:end-2,k))<0)+1;
            
            kp=find(CC_matSlave(ks,k)-CC_matSlave(ks+1,k)>0);
            ks=ks(kp);
            if length(ks)>nCallsXcorr-1
                
                [k1,k2]=sort(CC_matSlave(ks,k));
                ks=ks(k2(end-(nCallsXcorr-1):end));
                
                delay(:,k,jj)=ks;
                peak(:,k,jj)=k1(end-(nCallsXcorr-1):end);
                
                for j=1:nCallsXcorr
                    an=polyfit(-1:1,CC_matSlave(ks(j)-1:ks(j)+1,k)',2);
                    rdelay(j,k,jj)=ks(j)-an(2)/2/an(1);
                end;
            end
        end;
    end
    
    
    
    for jj=1:length(cc_matrix.slv)
        rdelay(:,:,jj)=(rdelay(:,:,jj)-(pha(jj)+1))/sample_freq*skip;
    end
    
    index=1;
    
    for(k=1:length(hyd(local_array(1)).detection.calls))
        hydro.sf2(k,:)=[hyd(local_array(1)).detection.calls(k).start_time,...
            hyd(local_array(1)).detection.calls(k).end_time];
    end;
    
    %%%%%%%  from here on, assume nm=1
    
    test0=squeeze(delay);
    k=find(test0);test0(k)=1;
    red_test=sum(test0,2);
    k=find(red_test>=min_pairs);  % don't consider time slices with few
    trim=k;
    test0=test0(k,:);  % delays than min required for a
    % localization
    
    k=find(trim+num_calls-1<=nOtherCalls);
    trim=trim(k);test0=test0(k,:); % avoid endpoint problem for rtime computation
    
    % start looping over ensembles with nontrivial entries
    
    for(jj0=1:length(trim))
        
        perccount(jj0,length(trim));
        jj=trim(jj0);
        
        if isfield(localize_struct.parm,'min_pairs')
            
            clear mt vec
            
            % Create columns (?) for cross correlations?
            rng=[length(local_array)-1:-1:localize_struct.parm.min_pairs];
            
            
            % For eeach slave hydrophone...
            for kl0=1:length(rng)
                
                kl=rng(kl0);
                
                % Create array with all combinations of slave MARUS (??)
                comb=nchoosek([2:length(local_array)],kl);
                
                [s1,s2]=size(comb);
                clear mt0 vec0
                
                for jl=1:s1
                    %   [jj0,comb(jl,:)-1]
                    
                    check=sum(test0(jj0,comb(jl,:)-1));
                    if check<kl % not enough delays
                        mt0(jl)=nan;
                        vec0(jl,:)=zeros(1,2)+nan;
                    else
                        
                        % there are enough delays so create the lsq grid
                        lsq=zeros(length(latgrid),length(longrid),kl);
                        
                        for(kk=1:kl)
                            lsq(:,:,kk)=(array_struct(array_number).toa_diff{comb(jl,kk)}...
                                -rdelay(1,jj,comb(jl,kk)-1)).^2;end
                        
                        test=sum(lsq,3);
                        [k1,k2]=min(test);[mt0(jl),kk2]=min(k1);
                        kk1=k2(kk2);
                        vec0(jl,:)=[latgrid(kk1),longrid(kk2)];
                        
                        % refine point unless near the edge of the grid
                        if kk1> 2 && kk1 < length(latgrid)-1
                            if kk2 > 2 && kk2 < length(longrid)-1
                                tma=test(kk1-2:kk1+2,kk2-2:kk2+2);
                                tmai=interp2(tma,5,'cubic');
                                [j1,j2]=min(tmai);[mt0(jl),jj2]=min(j1);
                                jj1=j2(jj2);
                                lat_loc=linspace(latgrid(kk1-2),latgrid(kk1+2),129);
                                lon_loc=linspace(longrid(kk2-2),longrid(kk2+2),129);
                                %                 mtest
                                %                 imagesc(lon_loc,lat_loc,tmai);pause
                                vec0(jl,:)=[lat_loc(jj1),lon_loc(jj2)];
                            end
                        end
                    end % test if delays make sense
                end % jl loop
                
                % best for the set just checked:
                [k1,k2]=min(mt0);
                mtest(kl0)=mt0(k2);vec(kl0,:)=vec0(k2,:);
            end % kl loop
            
            % now pick mtest as the best of these, allowing for natural
            % decrease of lsq error with fewer pairs (linear in #pairs?)
            
        else
            
            lsq=zeros(length(latgrid),length(longrid),length(local_array)-1);
            for(kk=1:length(local_array)-1)
                lsq(:,:,kk)=(array_struct(array_number).toa_diff{kk+1}-rdelay(1,jj,kk)).^2;end
            test=sum(lsq,3);
            [k1,k2]=min(test);[mtest,kk2]=min(k1);
            kk1=k2(kk2);
            vec=[latgrid(kk1),longrid(kk2)];
            
            % refine point unless near the edge of the grid
            if kk1> 2 && kk1 < length(latgrid)-1
                if kk2 > 2 && kk2 < length(longrid)-1
                    tma=test(kk1-2:kk1+2,kk2-2:kk2+2);
                    tmai=interp2(tma,5,'cubic');
                    [j1,j2]=min(tmai);[mtest,jj2]=min(j1);
                    jj1=j2(jj2);
                    lat_loc=linspace(latgrid(kk1-2),latgrid(kk1+2),129);
                    lon_loc=linspace(longrid(kk2-2),longrid(kk2+2),129);
                    %                 mtest
                    %                 imagesc(lon_loc,lat_loc,tmai);pause
                    vec=[lat_loc(jj1),lon_loc(jj2)];
                end
            end
            
        end % min_pairs test
        
        %      if mtest<cutoff  % A point satisfies cutoff criterion
        
        delays(index,:)=rdelay(1,jj,:);
        for k=1:length(local_array)-1;
            cross_score(index,k)=peak(1,jj,k);
        end
        times(index,:)=hydro.sf2(jj,:);
        rtimes(index)=mean(mean(hydro.sf2(jj:jj+num_calls-1,:)));
        dex(index)=jj; % remember when event happened
        
        if isfield(localize_struct.parm,'min_pairs')
            score(:,index)=mtest;
            coordinates(:,:,index)=vec;
        else
            score(index)=mtest;
            coordinates(:,index)=vec;end
        
        index=index+1;
        
        %      end % cutoff
        
        
    end % jj0 loop
    
    %
    localize_struct.hyd(local_array(1)).delays=delays;
    localize_struct.hyd(local_array(1)).cross_score=cross_score;
    localize_struct.hyd(local_array(1)).coord_time=times;
    localize_struct.hyd(local_array(1)).rtimes=rtimes;
    localize_struct.hyd(local_array(1)).dex=dex;
    localize_struct.hyd(local_array(1)).coordinates=coordinates;
    localize_struct.hyd(local_array(1)).score=score;
    
    
    
    
    
    
    
