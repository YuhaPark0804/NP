function [t, y] = odeRK4_student(odeFunc,a,b,h,y0)
% Description:
% % [TO-DO] Write Comment Section 
% % [TO-DO] Write Comment Section 
% % [TO-DO] Write Comment Section 
%
% Input variables:
% odeFunc:    Name of a function file that calculates dy/dt.  
% % [TO-DO] Write Comment Section 
% % [TO-DO] Write Comment Section 
% % [TO-DO] Write Comment Section 





    % Variable Initialization
    N = (b-a)/h;
    y=zeros(1,N+1);
    t=zeros(1,N+1);

    % Initial Condition
    t(1) = a;  y(1) = y0;
    

    % RK Design Parameters   
    alpha2 = 1/2;
    alpha3 = alpha2;
    beta21 = alpha2;
    beta32 = alpha2;

    alpha4 = 1;
    beta43 = alpha4;
    beta31 = 0;
    beta41 = beta31;
    beta42 = beta31;

    C1=1/6;
    C4=C1;
    C2=2/6;
    C3=C2;

    % ODE Solver
    for i = 1:N
        t(i+1) = t(i) + h;

        % [First-point Gradient] --------------------------------
        K1 = odeFunc(t(i),y(i));
 
        t2 = t(i) + alpha2 * h;
        y2 = y(i)+ beta21*K1*h;
        K2=odeFunc(t2,y2);
      
        % [Second-point Gradient]---------------------------
           
        t3 = t(i) + alpha3 * h;
        y3=y(i)+ beta31*K1*h +beta32*K2*h;
        K3=odeFunc(t3,y3);

        % [Second-point Gradient]---------------------------
           
        t4 = t(i) + alpha4 * h;
        y4 = y(i)+ beta41*K1*h + beta42*K2*h + beta43*K3*h;
        K4=odeFunc(t4,y4);

        % Estimate: y(i+1) using RK2  
        y(i+1)=y(i)+ (C1*K1 + C2*K2+ C3*K3+ C4*K4)*h;

        
    end

end  % END OF FUNCTION