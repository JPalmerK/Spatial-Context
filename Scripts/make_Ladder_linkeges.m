function chain = make_Ladder_linkeges(Corr_coef_map, localize_struct, ...
    master, fs, cutoff, time_cut)
% Function to create call groups on the correlation plots (from Glen)
% Inputs
%   fs - sample frequency 
%   cutoff - correlation probability coefficient
%   time_cut - maximum delay beteen calls to correlate

time = localize_struct.hyd(master).rtimes;


% Elapsed time since first call
time=(time-time(1))/fs;  % seconds now yes?
time0=time; % keep original copy for later indexing

[s1,s2]=size(Corr_coef_map);
nc=1; % cluster counter
clear chain dex ss

while s1>1
    flag=1;
    
    % call element index, start at 1 and examine columns
    elt=1; % first element of cluster is always upper left hand element
    index=[];
    cluster=[];
    
    
    % Find all calls that are sufficiently correlated with this one.
    while flag == 1
        % Add call to current cluster, retaining the time and index
        % of the call element
        cluster = [cluster, time(elt)]; 
        index   = [index,elt];
        
        % entire column for each element
        vec=Corr_coef_map(:,elt);
        
        % Only look at the elements below
        vec=vec(elt+1:end);
        
        % Look for correlations that exceed the threshold
        k=find(vec>cutoff);
        
        if length(k)==0
            flag=0;  % Nobody met the criterion
        else
            elt=k(1)+elt;  % Set next call to the first that meets criterion
        end
    end % flag test
    
    % now cut cluster at first large jump:
    k=find(diff(cluster)>time_cut);
    if length(k)>0
        k=k(1);
        % Remove everything after the first large gap
        cluster=cluster(1:k);
        index=index(1:k);
    end
    
    % aggregate clusters
    chain(nc).links=cluster;
    chain(nc).n = length(cluster);
    for j=1:length(cluster)
        chain(nc).index(j)=find(time0==cluster(j));
    end
    nc=nc+1;
    
    % filter out this cluster to reduce matrix
    [s1,s2]=size(Corr_coef_map);
    rs1=setdiff([1:s1],index);
    Corr_coef_map=Corr_coef_map(rs1,:);
    Corr_coef_map=Corr_coef_map(:,rs1);
    time=time(rs1);
    [s1,s2]=size(Corr_coef_map);
    
end % while s1 loop

chain = struct(chain);

end
