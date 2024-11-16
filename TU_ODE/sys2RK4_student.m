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
    alpha2 = 0.5;
    alpha3 = 0.5;
    beta21 = 0.5;
    beta32 = 0.5;

    alpha4 = 1;
    beta43 = alpha4;
    beta31 = 0;
    beta41 = beta31;
    beta42 = beta31;

    C1 = 1/6;
    C4 = C1;
    C2 = 2/6;
    C3 = C2;

    for i = 1:N
        t(i+1) = t(i) + h;
        
        % First-point Gradient
        dYdt = odeFunc(t(i),[y(i), v(i)]);
        Ky1 = dYdt(1);
        Kv1 = dYdt(2);
                
        % Second-point Gradient  
        t2 = t(i) + alpha2*h;
        y2= y(i) + beta21*Ky1*h;
        v2=v(i)+ beta21*Kv1*h;
    
        dYdt2=odeFunc(t2,[y2, v2]);
 
        Ky2 = dYdt2(1);
        Kv2 = dYdt2(2);

        % Third-point Gradient  
        t3 = t(i) + alpha3*h;
        y3 = y(i) + beta31*Ky1*h + beta32*Ky2*h;
        v3 = v(i) + beta31*Kv1*h + beta32*Kv2*h;
    
        dYdt3 = odeFunc(t3,[y3, v3]);
 
        Ky3 = dYdt3(1);
        Kv3 = dYdt3(2);

        % Forth-point Gradient  
        t4 = t(i) + alpha4*h;
        y4 = y(i) + beta41*Ky1*h + beta42*Ky2*h + beta43*Ky3*h;
        v4 = v(i) + beta41*Kv1*h + beta42*Kv2*h + beta43*Kv3*h;
    
        dYdt4 = odeFunc(t4,[y4, v4]);
 
        Ky4 = dYdt4(1);
        Kv4 = dYdt4(2);

        %out
        y(i+1) = y(i) + (C1*Ky1 + C2*Ky2 + C3*Ky3 + C4*Ky4)*h;  
        v(i+1) = v(i) + (C1*Kv1 + C2*Kv2 + C3*Kv3 + C4*Kv4)*h;
        
    end

end  % END OF FUNCTION