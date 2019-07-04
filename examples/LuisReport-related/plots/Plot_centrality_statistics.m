% Initialise arrays storing the ploted results
means_aggregated_hmm = zeros(3,3,3);                 % (Repetitions x Conditions x States)
means_aggregated_hsmm = zeros(3,3,3);

modes_aggregated_hmm = zeros(3,3,3);
modes_aggregated_hsmm = zeros(3,3,3);

means_states_hmm = zeros(27,5);
means_states_hsmm = zeros(27,5);

% Get statistic for each condition and repetition (9 conditions * 3 repetitions, concatenated as [9 9 9]).
for i = 1:27      
    
    % Keep track of which repetition are we in.
    if i < 10
        rep = 1;
    elseif i < 19
        rep = 2;
    else
        rep = 3;
    end
    
    % Get true-to-inferred states (Boolean) matching array.
    matching_hmm = results(i).HMM.true_to_inferred_mapping;
    matching_hsmm = results(i).HSMM.true_to_inferred_mapping;
    
    % Get the empirical distributions (reference/true and from HMM) and HSMM-modelled distribution.
    reference = get_empirical_lifetimes(results(i).hsmmtrue.states_seq, results(i).hsmmtrue.K, results(i).hsmmtrue.M);
    d_hmm = results(i).HMM.pd;
    d_hsmm = results(i).HSMM.model.sojourns;
    
    % Get reference means and modes from the reference (for each state).
    reference_means = get_mean_lifetime(reference);
    [~, reference_modes] = max(reference);
    
    % For each state, get mode and mean of the inferred distributions and calculate residuals with respect to the reference.
    for s = 0:4
        
        % Which inferred states are modelling the true state 's'?
        ms_hmm = find(matching_hmm(s+1,:));
        ms_hsmm = find(matching_hsmm(s+1,:));
        
        %%% For the HMM
        
        % Choose longer one if redundant states exist
        if length(ms_hmm) > 1 
            [~, modes] =  max(d_hmm(:,ms_hmm));
            [~, ind] = max(modes);
            ms = ms_hmm(ind);
        else
            ms = ms_hmm;
        end
        
        % Get means and modes
        mean_state = get_mean_lifetime(d_hmm(:,ms));
        [~, mode_state] = max(d_hmm(:,ms));

        means_states_hmm(i,s+1) = mean_state;
        
        % Add residuals to the aggregated residuals sum over states.
        means_aggregated_hmm(i) = means_aggregated_hmm(i) + abs(mean_state - reference_means(s+1));
        modes_aggregated_hmm(i) = modes_aggregated_hmm(i) + abs(mode_state - reference_modes(s+1));
        
        
        %%% Same for the HSMM
        
        % Choose longer one if redundant states exist
        if length(ms_hsmm) > 1 
            [~, modes] =  max(d_hsmm(:,ms_hsmm));
            [~, ind] = max(modes);
            ms = ms_hsmm(ind);
        else
            ms = ms_hsmm;
        end
        
        % Get means and modes
        mean_state = get_mean_lifetime(d_hsmm(:,ms));
        [~, mode_state] = max(d_hsmm(:,ms));

        means_states_hsmm(i,s+1) = mean_state;
        
        % Add residuals to the aggregated residuals sum over states.
        means_aggregated_hsmm(i) = means_aggregated_hsmm(i) + abs(mean_state - reference_means(s+1));
        modes_aggregated_hsmm(i) = modes_aggregated_hsmm(i) + abs(mode_state - reference_modes(s+1));
        
    end
end 

% Statistical testing comparing the two groups (H0: no differences).
%corrected_significance_means = compare_groups(means_aggregated_hmm, means_aggregated_hsmm)
%corrected_significance_modes = compare_groups(modes_aggregated_hmm, modes_aggregated_hsmm)

%disp(p_val)

average_means_states_hmm = sum(reshape(means_states_hmm, 3,3,3,5),3) / 3;
average_means_states_hsmm = sum(reshape(means_states_hsmm, 3,3,3,5),3) / 3;
figure
for i = 1:5
   subplot(2,3,i)
   heatmap(average_means_states_hmm(:,:,:,i));
   colorbar('off')
   title(strcat('State ', num2str(i)))
end
sgtitle('State-specific Mean Residuals, HMM')

figure
for i = 1:5
   subplot(2,3,i)
   heatmap(average_means_states_hsmm(:,:,:,i));
   colorbar('off')
   title(strcat('State ', num2str(i)))
end
sgtitle('State-specific Mean Residuals, HSMM')


% Normalize by K * repetitions (each point in the matrix is the sum of residuals of K states for each of the 3 repetitions).
means_aggregated_hmm = reshape(sum(means_aggregated_hmm,3) / (3*5), 3,3);
modes_aggregated_hmm = reshape(sum(modes_aggregated_hmm,3) / (3*5), 3,3);
    
means_aggregated_hsmm = reshape(sum(means_aggregated_hsmm,3) / (3*5), 3,3);
modes_aggregated_hsmm = reshape(sum(modes_aggregated_hsmm,3) / (3*5), 3,3);


figure
% Plot aggregated JSdiv
    % Axis
    STN = [2 5 10];
    Coupling = [0.3 0.5 0.8];
    Z_hmm_means = means_aggregated_hmm;
    Z_hsmm_means = means_aggregated_hsmm;
    
    Z_hmm_modes = modes_aggregated_hmm;
    Z_hsmm_modes = modes_aggregated_hsmm;
    
    
    %%% HMM
    
    % Means
    subplot(2,2,1)
    heatmap(STN, Coupling, Z_hmm_means);
    caxis([0,26])
    colorbar('off');
    title('HMM (means)')
    xlabel('STN')
    ylabel('Coupling')
    
    % Modes
    subplot(2,2,2)
    heatmap(STN, Coupling, Z_hmm_modes);
    caxis([0,26])
    colorbar('off');
    title('HMM (modes)')
    xlabel('STN')
    ylabel('Coupling')
    
    %%% HSMM
    
    % Means
    subplot(2,2,3)
    heatmap(STN, Coupling, Z_hsmm_means);
    caxis([0,26])
    colorbar('off');
    title('HSMM (means)')
    xlabel('STN')
    ylabel('Coupling')
    
    % Modes
    subplot(2,2,4)
    heatmap(STN, Coupling, Z_hsmm_modes);
    caxis([0,26])
    colorbar('off');
    title('HSMM (modes)')
    xlabel('STN')
    ylabel('Coupling')

    sgtitle('Mean residuals of lifetime centrality estimates (aggregated over states)')

