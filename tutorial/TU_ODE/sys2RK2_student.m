function [t, y, v] = sys2RK2_student(odeFunc,a,b,h,yINI,vINI)
% It solves a system of two first-order initial value ODEs using
% second-order Runge-Kutta method.
% The independent variable is x, and the dependent variables are y and z.
%
% Input variables:
% odeFunc   Name of a function file that calculates [dy/dt ; dv/dt].
% a      The first value of x.
% b      The last value of x.
% h      The size of a increment.
% yINI     The initial value of y.
% vINI     The initial value of v.
%
% Output variable:
% t      A vector with the t coordinate of the solution points.
% y      A vector with the y coordinate of the solution points.
% v      A vector with the v coordinate of the solution points.

% NOTE:
% [dYdt] = odeFunc(t,vecY)
%       Input:  vecY=[y ; z] 
%       Output: dYdt=[dydt ; dzdt] 

   % Variable Initialization
    N = (b-a)/h;    
    t=zeros(1,N+1);    
    y=zeros(1,N+1);
    v=zeros(1,N+1);

    % Initial Condition
    y(1) = yINI;
    v(1) = vINI;
    t(1)=a;

    % RK Design Parameters    
    alpha=1;
    beta=alpha;
    C1=0.5;
    C2=1-C1;

    for i = 1:N
        t(i+1) = t(i) + h;
        
        % First-point Gradient
        dYdt = odeFunc(t(i),[y(i), v(i)]);
        Ky1 = dYdt(1);
        Kv1 = dYdt(2);
                
        % [Second-point Gradient]  
        t2 = t(i) + alpha*h;
        y2= y(i) + beta*Ky1*h;
        v2=v(i)+ beta*Kv1*h;
        % Second-point Gradient
        dYdt2=odeFunc(t2,[y2, v2]);
 
        Ky2 = dYdt2(1);
        Kv2 = dYdt2(2);
          
        y(i+1) = y(i) + (C1*Ky1 + C2*Ky2)*h;  
        v(i+1) = v(i) + (C1*Kv1 + C2*Kv2)*h;
        
    end

end  % END OF FUNCTION