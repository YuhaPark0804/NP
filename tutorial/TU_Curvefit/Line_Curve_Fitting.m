clc;clear all; close all;

T=[30	40	50	60	70	80];
P=[1.05 1.07    1.09	1.14	1.17	1.21];

% MATLAB function of curvefitting
% z=[a1, a0]; 
Z=polyfit(T,P,1)

% yopt=Z(1).*x+Z(2);
x=30:10:150;
yopt=polyval(Z,x); 