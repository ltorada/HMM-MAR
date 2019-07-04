function P = matrix_multip_logspace(A, B)
    
    % Matrix multiplication in logspace.
    
    size_A = size(A);
    size_B = size(B);
    
    P = zeros(size_A(1), size_B(2));
    
    for j = 1:size_B(2)
        P(:,j) = nl_logsumexp(A + B(:,j)', 2);
    end
    
end