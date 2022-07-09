function [CE,SR,SD, mu] = Performance_Measure(Portfolio_return)
%PERFORMANCE_MEASURE Summary of this function goes here
%{
Input - 
  Portfolio_return - column vector of returns 

Output
CE - certainity equivalent - 1x1
SR - Sharpe Ratio - 1x1
SD - standard deviation - 1x1



%}

gamma = 2; % from lecture slides from Prof Pohlmeier

mu = mean(Portfolio_return);
SD = std(Portfolio_return);

CE = mu - (gamma/2)* (SD^2);
SR = mu/SD;







end

