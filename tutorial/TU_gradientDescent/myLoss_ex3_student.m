function E = myLoss_ex3_student(W,X,Y)

    %% Initialization
    a0=W(1);
    a1=W(2);
    m=length(X);
    E=0;
    dE=0;
    % sum_r = 0;

    %% Loss Function
    % Gradident of Loss Function 
    % E = (1/m) * sum(y-(a1*X+a0))^2;

    for i = 1:m
        %r = (Y(i)-(a1*X(i)+a0));
        %dE=dE+r^2;
        dE = dE + (Y(i)-(a1*X(i)+a0))^2;
    end

    % E = dE/m;
    E = (1/m) * dE;

end



    