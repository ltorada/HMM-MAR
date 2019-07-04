function m = get_mean_lifetime(lifetime_distribution)

    % Gets the means of state lifetime distributions.
    %
    % INPUT
    % -----
    %
    % lifetime_distribution     (M x S) matrix containing a lifetime distribution per column.
    %
    %
    % OUTPUT
    % ------
    %
    %  m                        (1 x S) vector containing the state-specific means.
    %
    %
    % Author: Luis Torada Aguilella
    
    
    M = size(lifetime_distribution, 1);
    S = size(lifetime_distribution, 2);
    
    m = sum(repmat(1:100, S, 1)' .* lifetime_distribution, 1) ./ sum(lifetime_distribution,1);
    
end

