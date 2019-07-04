function pd = get_empirical_lifetimes(vpath, S, M)

% Get empirical distribution from Viterbi path (S = number of states, M =
% maximal lifetime duration allowed).
%
%  
% Author: Luis Torada Aguilella
    

    for j = 1:S
        if ismember(j, vpath)
                
            % Trick to get the lifetimes fast
            vpath_j = (vpath==j);   
            d = [true, diff(vpath_j) ~= 0, true];
            lifetimes = diff(find(d));
            if vpath(1) == j
                lifetimes = lifetimes(1:2:end);
            else
                lifetimes = lifetimes(2:2:end);
            end
                
            pd(:,j) = pdf(fitdist(lifetimes','Kernel','BandWidth',1), 1:M);
        end
    end
    
end   