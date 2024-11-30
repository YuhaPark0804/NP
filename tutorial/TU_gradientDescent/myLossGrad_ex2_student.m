function dL= myLossGrad_ex2_student(W)

    x=W(1);
    y=W(2);
    dL= [0;0];

    % [TO-DO] Gradident of Loss Function 
    dL= [6*(x-2);
        2*(y-2)];

end
