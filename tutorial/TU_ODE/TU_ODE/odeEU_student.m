function [t, y] = odeEU_student(odeFunc,a,b,h,y0)
% odeED solves a first order initial value ODE using Euler Explicit Method
%
%
% Input variables:
% odeFunc    Name of a function file that calculates dy/dt.  
% a      The first value of t.
% b      The last value of t.
% h      Step size.
% yIni     The value of the solution y at the first point (initial value).
% Output variable:
% t      A vector with the t coordinate of the solution points.
% y      A vector with the y coordinate of the solution points. 


    % Variable Initialization
    N = (b-a)/h;
    y=zeros(1,N+1);
    t=zeros(1,N+1);
    
    % Initial Condition
    t(1) = a;  
    y(1) = y0;
    
    
    % ODE Solver
    for i = 1:N
        t(i+1) = t(i) + h;               
        y(i+1) = y(i) + odeFunc(t(i),y(i))*h;
    end


end  % END OF FUNCTION