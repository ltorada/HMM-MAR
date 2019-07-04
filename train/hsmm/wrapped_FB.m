function [Gamma, Xi, hmm] = wrapped_FB(hmm, Likelihood)
    
    %
    % DESCRIPTION:
    % -----------
    % FB replacing the nodescluster() function in hmminference.m .
    % nLL can be obtained by calling obslike() and then -log(ans). 
    %
    %
    %
    % INPUTS:
    % ------
    % hmm                   Structure with the initialized model (note: in addition to the requirements
    %                         for the hmm version, it must also contain: hmm.M, hmm.sojourns, ).
    % nLL                   -log likelihood of the data for each state (size = T x S).
    %
    %
    % OUTPUTS:
    % -------
    % Gamma                 Posterior distribution over states at each time point (size = T x S).
    % Xi                    Transition matrix at each time point given the data (size = T x S x S).                          
    % hmm                   hmm with the updated hmm.sojourns .
    %
    %
    % Author: Luis Torada Aguilella.
    %
    
    
 %% Extract initial conditions from the hmm structure
 
    %disp("### The HSMM FB is being used ###")
    
    S = hmm.K;                                          % Number of states.
    T = length(Likelihood);                                      % Number of time points.
    
    if isfield(hmm.train, "M")
        hmm.M = hmm.train.M;
    else
        hmm.M = 100;
    end
    M = hmm.M;                                          % Maximal lifetime allowed.
    
    if isfield(hmm.train, "sojourns_init")              % Sojourn distributions. (Must be (T x S)).
        hmm.sojourns = hmm.train.sojourns_init;
    elseif ~isfield(hmm, "sojourns")
        hmm.sojourns = 1/M * ones(M,hmm.K);             % Initialize with uniform distribution by default.
    end
    d = hmm.sojourns; 
    
    alpha0 = hmm.Pi;                                    % Prior for alpha0.
    
    Aij = hmm.P;                                        % Transition matrix.
    % Make P zero-diagonal for the hsmm FB.
    distributed_diag_prob = diag(Aij) ./ (hmm.K-1);
    Aij = Aij + distributed_diag_prob;
    Aij = Aij - diag(Aij) .* eye(hmm.K);
    
    hmm.P = Aij;
    
 %% Glue to C++
    
    tau = int64(T);
    censoring = int64(2);
    J = int64(S);
    piPara = alpha0;
    pPara = reshape(Aij', [1, S*S]); 
    dPara = reshape(d, [1,M*S]);  % p(d)
    pdfPara = reshape(Likelihood, [1, T*S]);       % likelihood
   
    
    % Initialize output variables
    F    = zeros(1, J * tau);
    L    = zeros(1, J * tau);
    G    = zeros(1, J * tau);
    L1   = zeros(1, J * tau);
    N    = zeros(1, tau);
    Norm = zeros(1, J * tau);
    eta  = zeros(1, J * M);
    xi   = zeros(1, J * M);
    err = int64(0);
    
    % Call MEX function 
    [~, L, ~, ~, ~, ~, ~, xi, ~] = FB(censoring, tau, J, M, dPara, pPara, piPara, pdfPara, F, L, G, L1, N, Norm, eta, xi, err);
    
    %assert(sum(L) ~= 0, 'L is a zeros vector.')
 %% Glue back to Matlab   
 
    % Reshape returned arrays
    gamma_tj = reshape(L, [tau, J]);         % Posterior over states (gamma).
    estimated_d = reshape(xi, [M, J]);       % Estimated sojourn distributions.
    
    norm = sum(estimated_d,1);
    if norm ~= 0
        estimated_d = estimated_d ./ norm;
    end
    hmm.sojourns_before_smoothing = estimated_d;
    [~, hmm.sojourn_modes_before_smoothing] = max(estimated_d);
    hmm.sojourn_means_before_smoothing = get_mean_lifetime(estimated_d);
    
 %% Smothing of p(d) and glue to Diego's code
 
    % Sample and smooth
 
        for j = 1:S
            %phat = gamfit(estimated_d(:,j));
            %hmm.sojourns(:,j) = pdf("Gamma", phat(1), phat(2), 1:M) ./ sum( pdf("Gamma", phat(1), phat(2), 1:M));
            if sum(estimated_d,1) ~= 0
                y = randsample(1:M,10000,true,estimated_d(1:M,j));
                pd = fitdist(y', 'Kernel', 'BandWidth',1);
                %pd = fitdist(y', 'Poisson');
                norm = sum(pdf(pd,1:M));
                %if norm == 0
                %hmm.sojourns(:,j) = pdf(pd,1:M);
                %hmm.sojourns(:,j) = 1/M * ones(M,hmm.K);
                %else
                hmm.sojourns(:,j) = pdf(pd,1:M) ./ norm;
                %hmm.sojourns(:,j) = pdf(pd,1:M);
                hmm.sojourn_means(j) = mean(pd);
                %end
            else
                hmm.sojourns(:,j) = 1/M * ones(M,1);
                hmm.sojourn_means(j) = M/2;
            end
        end
        
     %hmm.sojourns = estimated_d;

    
    Gamma = permute(gamma_tj', [2,1]);
    Gamma(Gamma < 0) = 0;
    Gamma = Gamma + realmin;
    
    % Eijt is the probability of a super-state ending at a given time and
    % transitioning to another super-state (super-state = sequence of equal
    % states (j,d) ), NOT the probability of a transition in general. So we
    % need to calculate Xi de novo, which is the joint probability of a parent
    % and child responsibilities, normalized with the data.
    
    Gamma_previous = Gamma(1:end-1, :);
    Gamma_next = Gamma(2:end, :);
    Xi = nan(T-1, S*S);
    for t = 1:T-1
        [i,j] = ndgrid(Gamma_previous(t,:), Gamma_next(t,:));
        all_transition_prob_combinations = [i(:),j(:)];
        Xi(t,:) = prod(all_transition_prob_combinations, 2)';
    end
    
    Xi(Xi<0) = 0;
    Xi = Xi + realmin;


end
