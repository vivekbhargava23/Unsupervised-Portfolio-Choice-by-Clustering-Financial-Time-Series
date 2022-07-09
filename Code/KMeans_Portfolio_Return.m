function [KMeans_Portfolio_return_OneByN_Model1,KMeans_Portfolio_return_OneByN_Model2, KMeans_Portfolio_Naive, KMeans_Portfolio_return_SR_Model3,TC_KMeans_Portfolio_return_OneByN_Model1,TC_KMeans_Portfolio_return_OneByN_Model2,TC_KMeans_Portfolio_return_SR_Model3] = KMeans_Portfolio_Return(Excess_Returns,T_obs,number_of_clusters)
%ANALYSIS Summary of this function goes here
%{  
INPUT - 
Excess_Return - TxN matrix with T observations for N assets
T_obs - T1x1 vector - usually 60x1 for 5 years analysis window


%}



[T,N] = size(Excess_Returns);



for i = 1:T-T_obs %############################## Check this value - T-T_obs
    
    X = Excess_Returns(i:T_obs+i-1,:); 
    
    %}
    clusters_final = kmeans(X',number_of_clusters,'Distance','sqeuclidean','MaxIter',100,'Replicates',5);
    
    Cluster_storage(:,i) = clusters_final;
    
    % K-means Clustering
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    weights_Model1 = OneByN_Model1(X,clusters_final);
    KMeans_Portfolio_return_OneByN_Model1(i,:) = weights_Model1 * Excess_Returns(T_obs+i,:)';
    
    KMeans_Portfolio_Naive(i,:) = ones(1,N)*(1/N)* Excess_Returns(T_obs+i,:)';
    
    weights_Model2 = OneByN_Model2(X,clusters_final);
    KMeans_Portfolio_return_OneByN_Model2 (i,:) = weights_Model2 * Excess_Returns(T_obs+i,:)';
    
    
    weights_SR_Model3 = SR_Model3(X,clusters_final);
    KMeans_Portfolio_return_SR_Model3(i,:) = weights_SR_Model3 * Excess_Returns(T_obs+i,:)';
    
    
    
    
    % Transaction costs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i==1
        c=0;
    else
        c= 0.0005; % Transaction cost 
    end
    
    % Model 1
    Store_weights_Model1(i+1,:) = weights_Model1;
    TC_Model1 = c * (sum(abs(Store_weights_Model1(i+1,:)-Store_weights_Model1(i,:))));
    TC_KMeans_Portfolio_return_OneByN_Model1(i,:) = KMeans_Portfolio_return_OneByN_Model1(i,:) - TC_Model1;
    
    
    % Model 2
    Store_weights_Model2(i+1,:) = weights_Model2;
    TC_Model2 = c * (sum(abs(Store_weights_Model2(i+1,:)-Store_weights_Model2(i,:))));
    TC_KMeans_Portfolio_return_OneByN_Model2(i,:) = KMeans_Portfolio_return_OneByN_Model2(i,:) - TC_Model1;
    
    
    % Model Sharpe ratio
    Store_weights_SR_Model3(i+1,:) = weights_SR_Model3;
    TC_Model3 = c * (sum(abs(Store_weights_SR_Model3(i+1,:)-Store_weights_SR_Model3(i,:))));
    TC_KMeans_Portfolio_return_SR_Model3(i,:) = KMeans_Portfolio_return_SR_Model3(i,:) - TC_Model3;
    
    

end







end

