%%% For each noise level:

% Signal-to-noise (STN) definition: activation noise - background noise, given a covariance of 3/4 that difference.
%
% S = 5, with d means = 5, 10, 20, 30, 40.
%
% Noises: 10, 5, 3, 2. 
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%    GENERATE DATA     %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = struct();
results.HMM = struct;
results.HSMM = struct;
i = 1;

timepoints = 10000;

for STN = [2 5 10  2 5 10  2 5 10]
    for Coupling = [0.3 0.5 0.8]
    
        % Define true model.
        sim_model_definition;
        
        % Store conditions
        results(i).STN = STN;
        results(i).Coupling = Coupling;
        results(i).timepoints = timepoints;
        
        % Simualate data.
        [X_hsmm,Gammatrue_hsmm, states_seq] = simhsmm(hsmmtrue);
        results(i).hsmmtrue = hsmmtrue;
        results(i).hsmmtrue.states_seq = states_seq';
        
        % Fit HMM
        [results(i).HMM.model, results(i).HMM.gamma, ~, results(i).HMM.vpath] = run_hmmar(X_hsmm, hsmmtrue, "hmm", 0);
        
        % Fit HSMM.
        [results(i).HSMM.model, results(i).HSMM.gamma, ~, results(i).HSMM.vpath] = run_hmmar(X_hsmm, hsmmtrue, "hsmm", 0);
        %results(i).HSMM.vpath = results(i).HSMM.vpath + 1;
        
        % Get empirical distributions
        results(i).HMM.pd = get_empirical_lifetimes(results(i).HMM.vpath', hsmmtrue.K, hsmmtrue.M);
        results(i).HSMM.pd = get_empirical_lifetimes(results(i).HSMM.vpath, hsmmtrue.K, hsmmtrue.M);
        
        % Match states to the 'true' (simulated) states
        results(i).HMM.true_to_inferred_mapping = get_states_correspondence(results(i).HMM, results(i).hsmmtrue);
        results(i).HSMM.true_to_inferred_mapping = get_states_correspondence(results(i).HSMM, results(i).hsmmtrue);
        
        i = i+1;
    end
end

%save 'Data_240619_banwidth1_widerD.mat' results
%save 'Data_260619_geo.mat' results


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%    PLOTS     %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Generate figures
    Plot_lifetimes
    Plot_covariances
    Plot_accuracy
    Plot_JSdiv
    Plot_centrality_statistics
    Plot_states_matching
