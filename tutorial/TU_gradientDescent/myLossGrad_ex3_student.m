function dE= myLossGrad_ex3_student(W,X,Y)

    %% Initialization
    a0=W(1);
    a1=W(2);
    m=length(X);
    dE=[0;0];
    sum_r = 0;
    sum_xr = 0;

    %% Loss Function
    % Gradident of Loss Function 

    for i = 1:m
        %r = (Y(i)-(a1*X(i)+a0));
        %dE=dE+[r; X(i)*r];

        sum_r = sum_r + (Y(i)-(a1*X(i)+a0));
        sum_xr = sum_xr + X(i) * (Y(i)-(a1*X(i)+a0));
    end

    %dE = -2./m.*(dE);
    dE = [(-2/m)*sum_r;
         (-2/m)*sum_xr];
    
end