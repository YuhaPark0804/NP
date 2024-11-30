function J=myJacobEx1_student(Z)
    % z(1)=x; % z(2)=y
    x = Z(1);
    y = Z(2);
    J=zeros(2,2);

    % Define Jacobian
    J=[-1/4*(exp(x/2)-exp(-x/2)), 1;
        18*x,                 50*y];     % 2x2 matrix for this exercise

end