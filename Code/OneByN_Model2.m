function [weights_assets,Total_clusters] = OneByN_Model2(X,clusters)
%ONEBYN_MODEL2 Summary of this function goes here
%   Detailed explanation goes here
%{   
INPUTS - 
X - TxN matrix where T is diff time period and N is the number of assets
clusters - is Nx1 vector of cluster number assigned to each N assets


%}

[T1,N1] = size(X);

train_size = .8;
test_size = 1-train_size;

train_observations = ceil(train_size * T1);



% Calculating weights for each cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Total_clusters = size(unique(clusters),1);

% finding indexes for each cluster group

index_clusters = zeros(Total_clusters,N1);


for i = 1:Total_clusters
    
    idx = find(clusters==i);
    number_observation_cluster(i)=size(idx,1);
    index_clusters(i,idx) = idx;
    
    T_sub = train_observations;
    subcluster_portfolio_return = zeros(T1-T_sub,1);
    
    
    for j = 1:T1-T_sub
        
        X_train = X(j:T_sub+j-1,idx);
        
        % Calculating GMVP weights for the subcluster
        weights_GMVP_subcluster = GMVP_weights(X_train);
        
        subcluster_portfolio_return(j,:) = weights_GMVP_subcluster' * X(T_sub+j,idx)';
        
        
        
    end
    % collating the returns on the out of sample data for each cluster
    cluster_return_testing_data(:,i) = subcluster_portfolio_return;
    Certainity_Equivalent_clusters(i) = Performance_Measure(subcluster_portfolio_return);
    
    
    
    
end




% Selecting the best clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cluster_threshold = 0.5; % select 50% of the clusters
a = ceil(Total_clusters * cluster_threshold);
Sorted_CE = sort(Certainity_Equivalent_clusters,'descend');
critical_value = Sorted_CE(a);
%critical_value = median(Certainity_Equivalent_clusters)

index_clusters_chosen = find(Certainity_Equivalent_clusters>=critical_value);


for ii = 1:size(index_clusters_chosen,2)
    
    
    index = find(clusters==index_clusters_chosen(ii));
    
    XX = X(:,index);
    Mean_X_selected_cluster(:,ii)  = mean(XX,2);
    
    
    
end

% Calculating the weights of the selected clusters using GMVP
w_GMVP = GMVP_weights(Mean_X_selected_cluster);

% Putting the cluster weights in the corresponding array
w_clusters = zeros(1,Total_clusters);
w_clusters(index_clusters_chosen) = w_GMVP;


% Calculating weights for each asset in individual cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

weights_assets = zeros(1,N1);

for i = 1:Total_clusters
    A = index_clusters(i,:);
    B = unique(A);
    C = setdiff(B,[0]);
    
    
    weights_assets(1,C) = w_clusters(i) * (1/number_observation_cluster(i));
    
end



end

