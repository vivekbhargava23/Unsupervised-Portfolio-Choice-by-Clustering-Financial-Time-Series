%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seminar - AEML
Submitted by -
 - Vivek Bhargava (01/1025699)

No. of functions used - 8

Details - 
1. Analysis
    a. Portfolio_Return - measure the portfolio returns for Hierarchcical
       clustering (with and without transaction costs)
    b. KMean_Portfolio_Return - measure the portfolio returns for K-means
       clustering (with and without transaction costs)

2. Asset Selection Models
    a. OneByN_Model1 - proposed model 1
    b. OneByN_Model2 - proposed model 2
    a. SR_Model3     - proposed model 3


3. Supplementary Function
    a. Performance_Measure - computes four merformance metric - CE, SR, std
       dev and Average Return
    b. GMVP_weights - computes the weights using the GMVP strategy
    c. Distance_Linkage - calculates the distance metric for hierarchical
       clustering

Other details -

  Data description
  TxN where T = 443 and N = 205 assets

Variations - 

- For changing distance criterion in hierarchical clustering change the
  values in Distance_Linkage.m file line 19 as 
  'spearman' or
  'euclidean'



  ------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
%
clear all;
close all;
clc
%

% Part1 - Reading Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts = detectImportOptions('sp500_price.csv');
opts = setvaropts(opts,'Date','InputFormat','MM/dd/uuuu');
dataset = readtable('sp500_price.csv',opts);
Date = dataset(:,1);

% Pre-Processing Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Table_data = dataset(:,~any(ismissing(dataset),1));
Table_data = table2array(Table_data(:,2:size(Table_data,2)-1));
 
[n,k] = size(Table_data);

% Calculating Returns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:n-1
    Return(i,:) = (Table_data(i+1,:) - Table_data(i,:))./ Table_data(i+1,:);
end

Risk_free_rate_percent = table2array(readtable('risk_free_rates.csv'));
Risk_free_rate = Risk_free_rate_percent(:,2) * 0.01;


% Calculating excess Return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Excess Return = return - monthly US treasury Rate
Excess_Return = Return - Risk_free_rate(2:size(Risk_free_rate,1),:);

% Excess_Return = Excess_Return(:,31:60);
[T,N] = size(Excess_Return);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Part 2 - Calculating Portfolio Performance for Hierarchical Clustering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_obs = 60; 


cluster_range = 10:2:68;
warning('off') 



for i = 1:size(cluster_range,2)
    i
    [Portfolio_return_OneByN_Model11,Portfolio_return_OneByN_Model22, Portfolio_Naive33,Portfolio_return_SR_Model33,TC_Portfolio_return_OneByN_Model1,TC_Portfolio_return_OneByN_Model2,TC_Portfolio_return_SR_Model3] = Portfolio_Return(Excess_Return,T_obs,cluster_range(i));
    
    %Calculating performance for one by n - Model 1
    Portfolio_return_OneByN_Model1(:,i) = Portfolio_return_OneByN_Model11;
    [CE1,SR1,SD1,mu1] = Performance_Measure(Portfolio_return_OneByN_Model11);
    Performance_Model1(:,i) = [CE1,SR1,SD1,mu1];
    
    %Calculating performance for one by n - Model 2
    Portfolio_return_OneByN_Model2(:,i) = Portfolio_return_OneByN_Model22;
    [CE2,SR2,SD2,mu2] = Performance_Measure(Portfolio_return_OneByN_Model22);
    Performance_Model2(:,i) = [CE2,SR2,SD2,mu2];
    
    %Calculating performance for Naive-Portfolio
    Portfolio_Naive(:,i) = Portfolio_Naive33;
    [CE_N,SR_N,SD_N,mu_N] = Performance_Measure(Portfolio_Naive33);
    Performance_Model_Naive(:,i) = [CE_N,SR_N,SD_N,mu_N]; 
    
    %Calculating performance for Cluster with Sharp Ratio Asset Selection
    Portfolio_return_SR_Model3(:,i) = Portfolio_return_SR_Model33;
    [CE_SR,SR_SR,SD_SR,mu_SR] = Performance_Measure(Portfolio_return_SR_Model33);
    Performance_Model_SR_cluster(:,i) = [CE_SR,SR_SR,SD_SR,mu_SR];
    
    
    
    % Calculating returns with Transaction Charges
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % One by n - model 1 with Transaction Costs
    Portfolio_return_OneByN_Model1_transaction(:,i) = TC_Portfolio_return_OneByN_Model1;
    [CE1_TC,SR1_TC,SD1_TC,mu1_TC] = Performance_Measure(TC_Portfolio_return_OneByN_Model1);
    Performance_Model1_transaction(:,i) = [CE1_TC,SR1_TC,SD1_TC,mu1_TC];
    
    
    % One by n - Model 2 with Transaction Costs
    Portfolio_return_OneByN_Model2_transaction(:,i) = TC_Portfolio_return_OneByN_Model2;
    [CE2_TC,SR2_TC,SD2_TC,mu2_TC] = Performance_Measure(TC_Portfolio_return_OneByN_Model2);
    Performance_Model2_transaction(:,i) = [CE2_TC,SR2_TC,SD2_TC,mu2_TC];
    
    
    % SR Model - with Transaction Costs
    Portfolio_return_SR_Model3_transaction(:,i) = TC_Portfolio_return_SR_Model3;
    [CE_SR_TC,SR_SR_TC,SD_SR_TC,mu_SR_TC] = Performance_Measure(TC_Portfolio_return_SR_Model3);
    Performance_Model_SR_cluster_transaction(:,i) = [CE_SR_TC,SR_SR_TC,SD_SR_TC,mu_SR_TC];
    
    
    
    

end


% Plotting the graphs
%{
Use
kk=1 : for plotting CE
kk=2 : for plotting Sharpe-ratio
kk=3 : for plotting Standard Deviation
kk=4 : for plotting Average Return
%}



% Calculating returns without Transaction Charges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(cluster_range,Performance_Model1(kk,:),'linewidth',2)
hold on
plot(cluster_range,Performance_Model2(kk,:),'linewidth',2)
hold on
plot(cluster_range,Performance_Model_Naive(kk,:),'linewidth',2)
hold on
plot(cluster_range,Performance_Model_SR_cluster(kk,:),'linewidth',2)
hold off
xline(30,'--r'); % change to 20 for spearman distance
xline(55,'--r');
a = get (gca, 'XTickLabel' );
set (gca, 'XTickLabel' , a, 'FontName' , 'Times' , 'fontsize' , 24)
legend({'Model1','Model2','Naive','Model Sharpe Ratio'})%,'Location','southeast')
title(' Hierarchical Model - excluding transaction cost (N=205)')
xlabel('Number of clusters')
ylabel('Average return')
ylim([-0.02 0.005])



T1 = array2table(Performance_Model1,'RowNames',{'Certainity Equivalent','Sharp Ratio','Std Dev','Average Return'})

T2 = array2table(Performance_Model2,'RowNames',{'Certainity Equivalent','Sharp Ratio','Std Dev','Average Return'})

T3 = array2table(Performance_Model_Naive,'RowNames',{'Certainity Equivalent','Sharp Ratio','Std Dev','Average Return'})

T4 = array2table(Performance_Model_SR_cluster,'RowNames',{'Certainity Equivalent','Sharp Ratio','Std Dev','Average Return'})





% Calculating returns with Transaction Charges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
h1 = plot(cluster_range,Performance_Model1_transaction(kk,:),'linewidth',2);
hold on
h2 = plot(cluster_range,Performance_Model2_transaction(kk,:),'linewidth',2);
hold on
h3 = plot(cluster_range,Performance_Model_Naive(kk,:),'linewidth',2);
hold on
h4 = plot(cluster_range,Performance_Model_SR_cluster_transaction(kk,:),'linewidth',2);
hold off
xline(30,':'); % change to 20 for spearman distance
xline(55,':');
legend([h1 h2 h3 h4 ],{'Model1','Model2','Naive Strategy','Model Sharpe Ratio'})%,'Location','southeast')
a = get (gca, 'XTickLabel' );
set (gca, 'XTickLabel' , a, 'FontName' , 'Times' , 'fontsize' , 24)
title(' Hierarchical Model including transaction cost(N=205)')
xlabel('Number of Clusters')
ylim([-0.02 0.005])


T1_transaction = array2table(Performance_Model1_transaction,'RowNames',{'Certainity Equivalent','Sharp Ratio','Std Dev','Average Return'})

T2_transaction = array2table(Performance_Model2_transaction,'RowNames',{'Certainity Equivalent','Sharp Ratio','Std Dev','Average Return'})

T3_transaction = array2table(Performance_Model_Naive,'RowNames',{'Certainity Equivalent','Sharp Ratio','Std Dev','Average Return'})

T4_transaction = array2table(Performance_Model_SR_cluster_transaction,'RowNames',{'Certainity Equivalent','Sharp Ratio','Std Dev','Average Return'})





% Part 3 - Calculating Portfolio Performance for K-means Clustering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%


for i = 1:size(cluster_range,2)
    i
    [Portfolio_return_OneByN_Model11,Portfolio_return_OneByN_Model22, Portfolio_Naive33,Portfolio_return_SR_Model33,TC_KMeans_Portfolio_return_OneByN_Model1,TC_KMeans_Portfolio_return_OneByN_Model2,TC_KMeans_Portfolio_return_SR_Model3] = KMeans_Portfolio_Return(Excess_Return,T_obs,cluster_range(i));
    
    %Calculating performance for one by n - Model 1
    Portfolio_return_OneByN_Model1(:,i) = Portfolio_return_OneByN_Model11;
    [CE1,SR1,SD1,mu1] = Performance_Measure(Portfolio_return_OneByN_Model11);
    KMeans_Performance_Model1(:,i) = [CE1,SR1,SD1,mu1];
    
    %Calculating performance for one by n - Model 1
    Portfolio_return_OneByN_Model2(:,i) = Portfolio_return_OneByN_Model22;
    [CE2,SR2,SD2,mu2] = Performance_Measure(Portfolio_return_OneByN_Model22);
    KMeans_Performance_Model2(:,i) = [CE2,SR2,SD2,mu2];
    
    %Calculating performance for Naive-Portfolio
    Portfolio_Naive(:,i) = Portfolio_Naive33;
    [CE_N,SR_N,SD_N,mu_N] = Performance_Measure(Portfolio_Naive33);
    KMeans_Performance_Model_Naive(:,i) = [CE_N,SR_N,SD_N,mu_N]; 
    
    %Calculating performance for Cluster with Sharp Ratio Asset Selection
    Portfolio_return_SR_Model3(:,i) = Portfolio_return_SR_Model33;
    [CE_SR,SR_SR,SD_SR,mu_SR] = Performance_Measure(Portfolio_return_SR_Model33);
    KMeans_Performance_Model_SR_cluster(:,i) = [CE_SR,SR_SR,SD_SR,mu_SR];
    
    
    
    
    % Calculating returns with Transaction Charges
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Kmeans_Portfolio_return_OneByN_Model1_transaction(:,i) = TC_KMeans_Portfolio_return_OneByN_Model1;
    [CE1_TC,SR1_TC,SD1_TC,mu1_TC] = Performance_Measure(TC_KMeans_Portfolio_return_OneByN_Model1);
    KMeans_Performance_Model1_transaction(:,i) = [CE1_TC,SR1_TC,SD1_TC,mu1_TC];
    
    
    % One by n - Model 2 with Transaction Costs
    Kmeans_Portfolio_return_OneByN_Model2_transaction(:,i) = TC_KMeans_Portfolio_return_OneByN_Model2;
    [CE2_TC,SR2_TC,SD2_TC,mu2_TC] = Performance_Measure(TC_KMeans_Portfolio_return_OneByN_Model2);
    KMeans_Performance_Model2_transaction(:,i) = [CE2_TC,SR2_TC,SD2_TC,mu2_TC];
    
    
    % SR Model - with Transaction Costs
    KMeans_Portfolio_return_SR_Model3_transaction(:,i) = TC_KMeans_Portfolio_return_SR_Model3;
    [CE_SR_TC,SR_SR_TC,SD_SR_TC,mu_SR_TC] = Performance_Measure(TC_KMeans_Portfolio_return_SR_Model3);
    KMeans_Performance_Model_SR_cluster_transaction(:,i) = [CE_SR_TC,SR_SR_TC,SD_SR_TC,mu_SR_TC];
    
    
    

end


% Plotting the graphs
%{
Use
kk=1 : for plotting CE
kk=2 : for plotting Sharpe-ratio
kk=3 : for plotting Standard Deviation
kk=4 : for plotting Average Return
%}



% Calculating returns without Transaction Charges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(cluster_range,KMeans_Performance_Model1(kk,:),'linewidth',2)
hold on
plot(cluster_range,KMeans_Performance_Model2(kk,:),'linewidth',2)
hold on
plot(cluster_range,KMeans_Performance_Model_Naive(kk,:),'linewidth',2)
hold on
plot(cluster_range,KMeans_Performance_Model_SR_cluster(kk,:),'linewidth',2)
hold off
legend({'Model1','Model2','Naive','Model Sharpe Ratio'})%,'Location','southeast')
a = get (gca, 'XTickLabel' );
set (gca, 'XTickLabel' , a, 'FontName' , 'Times' , 'fontsize' , 24)
title('Average return - KMeans Model excluding transaction cost (N=205)')
xlabel('Number of Clusters')
ylabel ('Average return')
ylim([-0.02 0.006])



T5 = array2table(KMeans_Performance_Model1,'RowNames',{'Certainity Equivalent','Sharp Ratio','Std Dev','Average Return'})

T6 = array2table(KMeans_Performance_Model2,'RowNames',{'Certainity Equivalent','Sharp Ratio','Std Dev','Average Return'})

T7 = array2table(KMeans_Performance_Model_Naive,'RowNames',{'Certainity Equivalent','Sharp Ratio','Std Dev','Average Return'})

T8 = array2table(KMeans_Performance_Model_SR_cluster,'RowNames',{'Certainity Equivalent','Sharp Ratio','Std Dev','Average Return'})




% Calculating returns with Transaction Charges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(KMeans_Performance_Model1_transaction(kk,:))
hold on
plot(KMeans_Performance_Model2_transaction(kk,:))
hold on
plot(KMeans_Performance_Model_Naive(kk,:))
hold on
plot(KMeans_Performance_Model_SR_cluster_transaction(kk,:))
hold off
legend({'Model1','Model2','Naive','Model Sharpe Ratio'})%,'Location','southeast')
a = get (gca, 'XTickLabel' );
set (gca, 'XTickLabel' , a, 'FontName' , 'Times' , 'fontsize' , 14)
title('Average Return Kmeans Model including transaction cost (N=205)')
ylim([-0.02 0.02])




T5_transaction = array2table(KMeans_Performance_Model1_transaction,'RowNames',{'Certainity Equivalent','Sharp Ratio','Std Dev','Average Return'})

T6_transaction = array2table(KMeans_Performance_Model2_transaction,'RowNames',{'Certainity Equivalent','Sharp Ratio','Std Dev','Average Return'})

T7_transaction = array2table(KMeans_Performance_Model_Naive,'RowNames',{'Certainity Equivalent','Sharp Ratio','Std Dev','Average Return'})

T8_transaction = array2table(KMeans_Performance_Model_SR_cluster_transaction,'RowNames',{'Certainity Equivalent','Sharp Ratio','Std Dev','Average Return'})






 
 
%}
%}
%}
