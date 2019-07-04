function nl_L = get_nl_L(nl_gamma_jt, T)
   
    % Get negative log-lokelihood of the data from the states posterior and
    % and the last timepoint. 
    
    nl_L = nl_logsumexp(nl_gamma_jt(:,T),1);
    
end