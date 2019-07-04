function S = matrix_add_logspace(A,B)
   
    % Amtrix addition in logspace
    
    size_A = size(A);
    size_B = size(B);

    % Account for cases where the number of columns differ by adding inf (zeros in Euclidean space).
    if size_B(1) < size_A(1)
        difference = size_A(1) - size_B(1);
        inf_patch = inf(difference, size_A(2));
        B = cat(1, B, inf_patch);
    end  
    
    % Do the addition one column at a time (assuming it is the larger dimension).
    S = nan(size_A);
    for col = 1:size_B(2)
        S(:,col) = nl_logsumexp([A(:,col), B(:,col)], 2);
    end
     
end