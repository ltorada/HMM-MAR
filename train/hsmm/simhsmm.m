function [X, Gammatrue, path] = simhsmm(hsmmtrue)

% Equivalent to simhmmmar.m for hsmm.
%
% INPUT
% -----
%
% hsmm              Structure containing the true model definition in the same format as the one outputed by hmmmar.
%
%
% OUTPUTS
% -------
%
% X                 (T x ndim) matrix containing the simulated observed data.
% path              Simulated sequence of hidden states (true Viterbi path).
% Gammatrue         Formatting of path as a posterior distribution with no uncertainty (used somewhere else. Size = (T x S).
%
%
% Author: Luis Torada Aguilella

    
    %% Markov Chain definition
    T = hsmmtrue.T;
    K = 1:hsmmtrue.K;                   % Number of states.    
    P = hsmmtrue.P;                     % Transition matrix.
    Pi = hsmmtrue.Pi;                   % Initial state probabilities.
    path = zeros(T,1);                  % Initialise sequence of states (length = T).
    Gammatrue = zeros(T,hsmmtrue.K);
    
    %% Simulate initial state and its permanence
    
    % Sample initial state from Pi
    sampled_state = randsample(K,1,true,Pi);
    
    % Initialise prev_n
    prev_n = 1;                                        % Index of the last simulated state.
    
    % Sample permanence from 1:T, update next_prev and include the state (1 + permanence) times in the record.
    if hsmmtrue.sojourns_form == "Geometric"
        permanence = round(random('Geometric', hsmmtrue.sojourns_parameters(sampled_state)));
    elseif hsmmtrue.sojourns_form == "Gamma"
        permanence = round(gamrnd(hsmmtrue.sojourns_parameters(sampled_state,1), hsmmtrue.sojourns_parameters(sampled_state,2), 1));
    end
    next_n = prev_n + permanence;                      % Index of the last repetition of the previously simulated state.
    if next_n > T
        next_n = T;
        path(prev_n:T) = sampled_state;            
    else
        path(prev_n:next_n) = sampled_state;           % Add repeated state until n (included).
    end
    
    next_state = sampled_state;
    %% Simulate the rest of states (up to N).
    
    while next_n < T
         
        %next_state = randsample(3,1,true,masked_weights);
        next_state = randsample(K,1,true,P(next_state,:));
        
        % Update prev_n.
        prev_n = next_n;
        
        % Sample permanence (as for the initial state), update next_n and include sequence of repetitions into the sequence of states record.
        if hsmmtrue.sojourns_form == "Geometric"
            permanence = round(random('Geometric', hsmmtrue.sojourns_parameters(next_state)));
        elseif hsmmtrue.sojourns_form == "Gamma"
            permanence = round(gamrnd(hsmmtrue.sojourns_parameters(next_state,1), hsmmtrue.sojourns_parameters(next_state,2), 1));
        end
        %permanence = round(random('Geometric', hsmmtrue.sojourns_parameters(next_state)));
        next_n = prev_n + permanence;
        if next_n>T
            next_n = T;                                          % Update current step
            path(prev_n:T) = next_state;                         % Add repeated state until the present
        else
            path(prev_n:next_n) = next_state;
        end
    end
    
    %% Simulate observations for the simulated sequence of states
    
    % Observation model definition: means of the Gaussians repressenting for state (keeping the variance constant)
    
    ndim = length(hsmmtrue.state(1).W.Mu_W);
    mu_dim_k = zeros(ndim, hsmmtrue.K);
    sigma_dim_k = zeros(ndim, ndim, hsmmtrue.K);
    for j = 1:hsmmtrue.K
        mu_dim_k(:,j) = hsmmtrue.state(j).W.Mu_W;
        %sigma_dim_k(:,j) = hsmmtrue.state(j).W.S_W;
        sigma_dim_k(:,:,j) = hsmmtrue.state(j).Cov;
    end
    
    % Initialise and generate sequence of observations.
    X = zeros(T, ndim);
    for state_index = 1:T
         X(state_index, :) = mvnrnd(mu_dim_k(:,path(state_index))', sigma_dim_k(:,:,path(state_index))', 1);
    end
  
    % Prepare Gammatrue
    for i = 1:T
      Gammatrue(i,path(i)) = 1; 
    end
end
    