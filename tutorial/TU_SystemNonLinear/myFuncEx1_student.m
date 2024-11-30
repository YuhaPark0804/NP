function F=myFuncEx1_student(Z)
    x = Z(1); 
    y = Z(2);

    F=zeros(2,1);

    F(1) = y-0.5*(exp(x/2)+exp(-x/2)) ;
    F(2) = 9*x^2+25*y^2-225;
end