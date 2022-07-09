function [weights_assets,Total_clusters] = OneByN_Model1(X,clusters)
%ONEBYN_MODEL1 Summary of this function goes here
%{   
INPUTS - 
X - TxN matrix where T is diff time period and N is the number of assets
clusters - is Nx1 vector of cluster number assigned to each N assets


%}
[T1,N1] = size(X);


% Calculating weights for each cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Total_clusters = size(unique(clusters),1);

% finding indexes for each cluster group

index_clusters = zeros(Total_clusters,N1);
Mean_return_cluster = zeros(T1,Total_clusters);


for i = 1:Total_clusters
    
    idx = find(clusters==i);
    number_observation_cluster(i)=size(idx,1);
    index_clusters(i,idx) = idx;
    
    XX = X(:,idx);
    Mean_return_cluster(:,i) = mean(XX,2);
end

w_GMVP_clusters = GMVP_weights(Mean_return_cluster);


% Calculating weights for each casset in individual cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

weights_assets = zeros(1,N1);

for i = 1:Total_clusters
    A = index_clusters(i,:);
    B = unique(A);
    C = setdiff(B,[0]);
    
    
    weights_assets(1,C) = w_GMVP_clusters(i) * (1/number_observation_cluster(i));
    
end


end

