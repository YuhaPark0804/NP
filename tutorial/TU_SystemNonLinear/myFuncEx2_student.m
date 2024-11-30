function F=myFuncEx2_student(Z)

    F=zeros(3,1);
    
    th=Z(1);
    dx=Z(2);
    dy=Z(3);
    
    % Initial Positions
    % P0
    P0x=0;          P0y=100; 
    % P1
    P1x=0;          P1y=-100;
    
    % Final Positions
    % P0
    P0x_new=50;     P0y_new=186.6025; 
    % P1
    P1x_new=150;    P1y_new=13.3975;
        
    % [TO-DO] Define Function F
    F = [cos(th)*P0x - sin(th)*P0y + dx - P0x_new;
         sin(th)*P0x + cos(th)*P0y + dy - P0y_new;
         cos(th)*P1x - sin(th)*P1y + dx - P1x_new];     % 3x1 vector for this exercise

end