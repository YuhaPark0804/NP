function J=myJacobEx2_student(Z)

J=zeros(3,3);
th=Z(1);
dx=Z(2);
dy=Z(3);
x0=0; y0=100; x1=0; y1=-100;

% [TO-DO] Define Jacobian
J=[(-sin(th))*x0 - (cos(th))*y0, 1, 0
   (cos(th))*x0 - (sin(th))*y0, 0, 1;
   (-sin(th))*x1 - (cos(th))*y1, 1, 0];     % 3x3 matrix for this exercise


end