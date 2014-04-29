function tau = compute_focusing_delays_linear_one_way(F,ac,c0)
    % returns the one-way focusing delay for a linear array 
    % F  - focus depth             
    % ac - points on the aperture [vector]
    % c0 - speed of sound  
        
    F = F(ones(1,length(ac)),:).';
    R = sqrt(ac.^2 + F.^2);
    tau = (R - F)/c0;
    
