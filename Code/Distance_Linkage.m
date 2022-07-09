function [Z_opt,max_cophen,char] = Distance_Linkage(X)
%DISTANCE&LINKAGE Summary of this function goes here
%   Detailed explanation goes here
%{
X = TxN matrix 
T - time series
N - no of assets

Output
Z_opt - Returns the linkage value which is a 3 column matrix
Column 1 & Column 2 shows the assets that are combined together
Column 3  - contains the linkage distance between the two clusters

%}


max_cophen = 0;

D1 = pdist(X','euclidean');
Y1 = squareform(D1);

% Defining different Linkages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single Linkage
%{
Z1 = linkage(Y1,'single');
c1 = cophenet(Z1,D1)

if c1>max_cophen
    max_cophen = c1;
    Z_opt = Z1;
    char = 'single'
end
%}

% Complete Linkage
Z2 = linkage(Y1,'complete');
c2 = cophenet(Z2,D1);

if c2>max_cophen
    max_cophen = c2;
    Z_opt = Z2;
    char = 'complete';
end

% Ward Linkage
Z3 = linkage(Y1,'ward');
c3 = cophenet(Z3,D1);


if c3>max_cophen
    max_cophen = c3;
    Z_opt = Z3;
    char = 'ward';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









end

