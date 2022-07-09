function [weights_assets,Total_clusters] = SR_Model3(X,clusters)
%SR_MODEL3 Summary of this function goes here
%   Detailed explanation goes here

%{   
INPUTS - 
X - TxN matrix where T is diff time period and N is the number of assets
clusters - is Nx1 vector of cluster number assigned to each N assets

%}
[T1,N1] = size(X);

Total_clusters = size(unique(clusters),1);


for i = 1:Total_clusters
    
    idx = find(clusters==i);
    number_observation_cluster(i)=size(idx,1);
    index_clusters(i,idx) = idx;
    
    XX_cluster = X(:,idx);
    SR_temp = mean(XX_cluster) ./ std(XX_cluster);
    [maximum Ind] = max(abs(SR_temp));
    
    asset_index_cluster(i) = idx(Ind);
    
end

X_asset_Sharp_Ratio = X(:,asset_index_cluster);

weights_SR_assets = GMVP_weights(X_asset_Sharp_Ratio);



weights_assets = zeros(1,N1);

for i = 1:Total_clusters
    
    weights_assets(:,asset_index_cluster(i)) = weights_SR_assets(i);

end



end

