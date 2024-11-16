function [t, y] = odeRK2_student(odeFunc,a,b,h,y0)
% Description: odeED solves a second order initial value ODE using 2nd Order Runge-Kutta
%
% Input variables:
% odeFunc:    Name of a function file that calculates dy/dt.  
% a      The first value of t.
% b      The last value of t.
% h      Step size.
% y0     Initial Condition of y(1)
%
% Output variable:
% t      A vector with the t coordinate of the solution points.
% y      A vector with the y coordinate of the solution points. 





    % Variable Initialization
    N = (b-a)/h;
    y=zeros(1,N+1);
    t=zeros(1,N+1);

    % Initial Condition
    t(1) = a;  y(1) = y0;
    

    % RK Design Parameters    
    alpha=1;
    beta=alpha;
    C1=0.5;
    C2=1-C1;

    % ODE Solver
    for i = 1:N
        t(i+1) = t(i) + h;

        % [First-point Gradient]
        K1 = odeFunc(t(i),y(i));
      
        % [Second-point Gradient]
        % calculate t2=t(i)+alpha*h   
        t2 = t(i) + alpha * h;

        % calculate y2=y(i)+ beta*K1*h   
        y2=y(i)+ beta*K1*h;

        % Calculate: K2   using t2 and y2  
        K2=odeFunc(t2,y2);

        % Estimate: y(i+1) using RK2   
        y(i+1)=y(i)+ (C1*K1 + C2*K2)*h;

        
    end

end  % END OF FUNCTION