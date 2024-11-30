function a_raw= unnormalize(Z,x)
    maxX = max(x);
    minX = min(x);

    a0_est = Z(1) - Z(2)*maxX/(maxX-minX);
    a1_est = Z(2)/(maxX-minX);
    
    a_raw=[a0_est a1_est];
end