function [vpath] = wrapped_Viterbi(hmm, Likelihood)
    
    %
    % DESCRIPTION:
    % -----------
    % HSMM version of the Viterbi algorithm wrapped from the hsmm R package and glue
    % to OHBA HMM-MAR in Matlab.
    %
    %
    % INPUTS:
    % ------
    % hmm                   Structure with the initialized model.
    % Likelihood            Likelihood of the data for each state (size = T x S).
    %
    %
    % OUTPUTS:
    % -------
    % vpath                 Most probable state timecourse (Viterbi path).
    %
    %
    % Author: Luis Torada Aguilella.
    %
    
    
 %% Extract initial conditions from the hmm structure
 
 
    %disp("### The HSMM FB is being used ###")
    
    S = hmm.K;                                       % Number of states.
    T = length(Likelihood);                          % Number of time points.
    
    if isfield(hmm.train, "M")
        hmm.M = hmm.train.M;
    else
        hmm.M = 500;
    end
    M = hmm.M;                                       % Maximal lifetime allowed.
    
    if isfield(hmm.train, "sojourns_init")           % Lifetime distributions. Must be (T x S).
        hmm.sojourns = hmm.train.sojourns_init;
    elseif ~isfield(hmm, "sojourns")
        hmm.sojourns = 1/M * ones(M,hmm.K);
        plot(hmm.sojourns)
    end
    d = hmm.sojourns; 
    
    alpha0 = hmm.Pi;                                 % Prior for alpha0
    
    Aij = hmm.P;                                     % Transition matrix (make zero-diagonal).
    distributed_diag_prob = diag(Aij) ./ (hmm.K-1);
    Aij = Aij + distributed_diag_prob;
    Aij = Aij - diag(Aij) .* eye(hmm.K);
    hmm.P = Aij;
    
 %% Glue to C++
    
    tau = int64(T);
    J = int64(S);
    piPara = alpha0;
    pPara = reshape(Aij', [1, S*S]); 
    dPara = reshape(d, [1,M*S]);  % p(d)
    pdfPara = reshape(Likelihood, [1, T*S]);       % likelihood
    
    % Initialize output variables
    %vpath = int64(zeros(1, tau));
    
    % Call MEX function 
    [vpath] =  Viterbi(tau, J, M, dPara, pPara, piPara, pdfPara);
    
    vpath = vpath + 1; % Since it was using indexing starting at 0.

end
