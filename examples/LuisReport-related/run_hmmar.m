function [fitted_structure, gamma, xi, vpath] = run_hmmar(X, hsmmtrue, model, d_init)
    
    % Run hmmmar.m sharing setting some options by copying them from the hsmmtrue
    % specifications (e.g. M, T, and p(d) if initialized explicitly).
    %
    % INPUTS
    % ------
    %
    % hsmmtrue             Structure used for simulating the data or that at least contains the fields T, M and K.
    % model                String specifying which model do you want to train: either "hmm" or "hsmm".
    % d_init               Binary variable specifying whether you have manually initialized p(d) by 
    %                      including a field d in hsmmtrue. If 0, a uniform distribution is used by default.
    %
    %
    % OUTPUTS
    % -------
    %
    % fitted_structure     Structure containing the trained model parameterization.
    % gamma                Posterior distribution over states at each timepoint.
    % xi                   Transition probabilities at each timepoint (based on gamma).
    % vpath                Viterbi path.
    %
    %
    % Author: Luis Torada Aguilella
    
    
    % Specify options:
    options = struct();
    options.K = hsmmtrue.K;
    options.Fs = 1;
    options.covtype = 'full';
    options.order = 0;
    options.DirichletDiag = 1;
    options.zeromean = 1;
    options.verbose = 1;
    options.dropstates = 0;
    options.useParallel = 0;
    options.initrep = 5;
    options.repetitions = 1;
    
    if model == "hsmm"
        
        % New options
         options.hsmm = 1;                     % 0 by default.
         options.M = hsmmtrue.M;               % 100 by default if not initialized.
         if d_init == 1
             options.sojourns_init = hsmmtrue.d;
         end
         [fitted_structure, gamma, xi, vpath] = hmmmar(X, hsmmtrue.T, options);
    
    elseif model == "hmm"
        
        options.hsmm = 0;
        
        [fitted_structure, gamma, xi, vpath] = hmmmar(X, hsmmtrue.T, options);
    end

end

