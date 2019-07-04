% Simulate hidden semi-markovian chain during the initialization of the state timecourse posterior.
% It uses a uniform distribution between [1, M].


%% Simulate initial state and its permanence
    
    path = zeros(T,1);  
    M = options.M;
    
    % Sample initial state from Pi
    sampled_state = randsample(K,1,true,Pi);
    
    % Initialise prev_n
    prev_n = 1;                                        % Index of the last simulated state.
    
    % Sample permanence from 1:T, update next_prev and include the state (1 + permanence) times in the record.
    permanence = randi([1,M], 1);

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
        permanence = randi([1,M], 1);
        
        next_n = prev_n + permanence;
        if next_n>T
            next_n = T;                                          % Update current step
            path(prev_n:T) = next_state;                         % Add repeated state until the present
        else
            path(prev_n:next_n) = next_state;
        end
    end
    
   
   for i = 1:T
      Gamma(i,path(i)) = 1; 
    end