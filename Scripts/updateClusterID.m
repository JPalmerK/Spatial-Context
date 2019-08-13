function Cluster_id= updateClusterID(obj)
            
           
            
            % Simply create an array with the predicted clusterid
            Cluster_id = zeros(size(obj.arrivalArray,1),1)+length(obj.chains)+1;
            
            % Step through each chain and grab the cluster
            for ii=1:length(obj.chains)
                Cluster_id(obj.chains(ii).index) = ii;
                
            end
            
        end