function [Gamma, Xi, hmm] = maxM_newd_vectorized_FB_V2(hmm, nLL)
    
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
    
    
 %% Initial conditions
    
     % Initialize veriables for initial conditions (include a padding of #M Inf).
 
    S = hmm.K;                                          % Number of states.
    T = length(nLL);                                    % Number of time points.
    
    if isfield(hmm.train, 'M')
        hmm.M = hmm.train.M;
    else
        hmm.M = hmm.train.M;
    end
    M = hmm.M;                                          % Maximal lifetime allowed.
    
    if ~isfield(hmm.train, 'sojourns_init')                        % Sojourn distributions. (Must be (T x S)).
        hmm.sojourns = -log(1/M) * ones(M,hmm.K);    
        %nl_d_tj = (nlog_d - nlog_norm_d)'; % In case we want to initialize with a parametric distribution.
    else
        hmm.sojourns = hmm.train.sojourns_init;
    end
    nl_d_tj = hmm.sojourns; 
      
    % Add padding to nLL
    nLL = cat(1,zeros(M,S), nLL, zeros(M,S));
    
    % Initial conditions for the Markov chain.
    nl_alpha0 = -log(hmm.Pi');                          % Prior for alpha0
    
    Aij = hmm.P;                                        % Transition matrix
    % Make P zero-diagonal for the hsmm FB.
    distributed_diag_prob = diag(Aij) ./ (hmm.K-1);
    Aij = Aij + distributed_diag_prob;
    Aij = Aij - diag(Aij) .* eye(hmm.K);
    nl_Aij = -log(Aij);
    
    %% Memory initialization
    
    nl_alpha_tj = inf(T+1+M,S);
    nl_alphaStar_tj = inf(T+1+M,S);
    
    nl_estimated_dj = inf(M,S);                % cummulative probability of being at (j,d).
  
    nl_beta_tj = inf(T+M,S);                   % Equivalent to forcing transitions after T to be of probability = 0.
    nl_betaStar_tj = inf(T,S);                 % probability of O_{t+1:T} given S_{[t+1} = j (that is, given that you transition to j).
    
    nl_gamma_tj = inf(T,S);
    gamma_tj = inf(T,S);   
 
    %% Get betas
    
    % Initialize last element
    nl_betaStar_tj(T,:) = -log(ones(1,S));                      % = not taking into account the evolution of the chain after t = T.
    nl_beta_tj(T,:) = -log(ones(1,S));
    
    step_b_tj = zeros(M,S);                                     % cummulative likelihood from t+step until t+M, for each step (step = index).
  
    % Recursively (backward in time) get the rest of the betas.
    for t = T-1:-1:1
        % bsxfun( @times, step_b_tj(1:end-1,:), nLL(M+1+t,:) )
        step_b_tj = cat(1,nLL(M+t+1,:), step_b_tj(1:end-1,:) + nLL(M+1+t,:));
        
        % for each possible transition from a given i to a given j, probability of being in j d timepoints and seeing the data.
        % corresponding to that interval of time.
        nl_betaStar_tj(t,:) = nl_logsumexp( nl_d_tj(1:M,:) + step_b_tj + nl_beta_tj(t+1:t+M,:), 1);    % size = (1 x S).
        
        % weight all possible transitions from each i to all j.
        nl_beta_tj(t,:) = matrix_multip_logspace( nl_Aij , nl_betaStar_tj(t,:)' );    
    end
    nl_beta_tj(T+1:end, :) = [];
    
    %% Get alphas
    
    % Initialize first element of alpha (t = 0).
    nl_alpha_tj(M+1,:) = matrix_multip_logspace(-log(hmm.Pi), nl_Aij);       
    total = matrix_multip_logspace(-log(hmm.Pi), nl_Aij);
    total(total == -inf) = inf;
    nl_alphaStar_tj(M+1,:) = total;
    
    % Initialize memory
    step_b_tj = zeros(M,S);                               % cummulative likelihood O_{t-M:t}   t-1 because the likelihood of t=0 doesn't count.  
    
    % Recursively (forwards in time) get the rest of the alphas.
    t = 1;
    step_b_tj(1,:) = nLL(M+t-1,:);
    for step = 2:M
        step_b_tj(step,:) = step_b_tj(step-1,:) + nLL(M+t-step,:);  % going from d=1 (index=1) until d=M.
    end
    
    for t = 2:T+1
        
        step_b_tj = cat(1, nLL(M+t-1,:) , step_b_tj(1:end-1,:) +  nLL(M+t-1,:));
        
        sum_space = nl_alphaStar_tj(M+t-1:-1:M+t-M,:) + nl_d_tj(1:M,:) + step_b_tj; % Careful: going from dj=t to dj=1
        
        % For updating p(d).
        nl_estimated_dj = matrix_add_logspace( nl_estimated_dj, sum_space + nl_beta_tj(t-1,:) );
        nl_estimated_dj(nl_estimated_dj == -inf) = inf;
        
        nl_alpha_tj(M+t,:) = nl_logsumexp(sum_space ,1);
        
        total = matrix_multip_logspace( nl_alpha_tj(M+t,:), nl_Aij );
        total(total == -inf) = inf;
        nl_alphaStar_tj(M+t,:) = total;
    
    end
    nl_alpha_tj(1:M+1,:) = [];
    
    %% Get Eijt
    nl_Eijt = inf(S,S,T-1);
    for t = 1:T-1
        nl_Eijt(:,:,t) = nl_alpha_tj(t,:)' + nl_Aij + nl_betaStar_tj(t,:);
    end
    
    %% Get gammas and back to Euclidean, normalized space
    
    nl_gamma_tj(end,:) = nl_alpha_tj(end,:)';
    
    nl_L = get_nl_L(nl_gamma_tj', T);
    
    gamma_tj(end,:) = exp(-(nl_gamma_tj(end,:) - nl_L));
    Eijt = exp(-(nl_Eijt - nl_L));
    
    estimated_dj = exp(-(nl_estimated_dj - nl_L));
    % Normalize
    estimated_dj = estimated_dj ./ sum(estimated_dj,1);
    
    
    % Sample and smooth
    for j = 1:S
        y = randsample(1:M,100000,true,estimated_dj(1:M,j));
        pd = fitdist(y', 'Kernel', 'BandWidth',4);
        hmm.sojourns(:,j) = -log(pdf(pd,1:M));
    end
    
    %hmm.sojourns = -log(estimated_dj);
    
    % Get the rest recursively
    for t = T-1:-1:1
        gamma_tj(t,:) = gamma_tj(t+1,:) + sum(Eijt(:,:,t) - Eijt(:,:,t)',2)';
    end 
    
    %% Adapt output format
  
    [m, inf_states_seq] = max(gamma_tj', [], 1);
    
    
    % In order to be merged with Diego's code:
    
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
