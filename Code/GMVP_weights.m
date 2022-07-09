function [w_GMVP] = GMVP_weights(returns)
%GMVP_WEIGHTS Summary of this function goes here
%   Detailed explanation goes here
%{
Inputs
returns - TxN Matrix
where T is the number of time periods
N is the number of assets
%}

[n k] = size(returns);
S = cov(returns);

vector_one=ones(size(S,1),1);

w_GMVP =(pinv(S)*vector_one)/(vector_one'*pinv(S)*vector_one);

end

