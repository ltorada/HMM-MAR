% Initialise arrays
JSdiv_aggregated_hmm = zeros(3,3,3); 
JSdiv_aggregated_hsmm = zeros(3,3,3);

JSdiv_states_hmm = zeros(27,5);
JSdiv_states_hsmm = zeros(27,5);

for i = 1:27      % For each simulated scenario
    
    
    matching_hmm = results(i).HMM.true_to_inferred_mapping;
    matching_hsmm = results(i).HSMM.true_to_inferred_mapping;
    
    % Keep track of which repetition are we in.
    if i < 10
        rep = 1;
    elseif i < 19
        rep = 2;
    else
        rep = 3;
    end
    
    % Get the empirical distributions for each condition.
    reference = get_empirical_lifetimes(results(i).hsmmtrue.states_seq, results(i).hsmmtrue.K, results(i).hsmmtrue.M);
    d_hmm = results(i).HMM.pd;
    d_hsmm = results(i).HSMM.model.sojourns;
    
    % For each state, get JS divergence.
    total_hmm = 0;
    total_hsmm = 0;
    for s = 1:5
        % Find state index from the HMM and HSMM matching the reference state s.
        ms_hmm = find(matching_hmm(s,:));
        ms_hsmm = find(matching_hsmm(s,:));
        
        % For the HMM
        if length(ms_hmm) > 1 % chose longer one
            [~, modes] =  max(d_hmm(:,ms_hmm));
            [~, ind] = max(modes);
            ms = ms_hmm(ind);
        else
            ms = ms_hmm;
        end
        JSdiv_states_hmm(i,s) = JSDiv(reference(:,s)', d_hmm(:,ms)');
        total_hmm = total_hmm + JSDiv(reference(:,s)', d_hmm(:,ms)');
        
        % For the HSMM
        if length(ms_hsmm) > 1 % chose longer one
            [~, modes] =  max(d_hsmm(:,ms_hsmm));
            [~, ind] = max(modes);
            ms = ms_hsmm(ind);
        else
            ms = ms_hsmm;
        end
        JSdiv_states_hsmm(i,s) = JSDiv(reference(:,s)', d_hsmm(:,ms)');
        total_hsmm = total_hsmm + JSDiv(reference(:,s)', d_hsmm(:,ms)');
    end

    % Now average-out the JS divergence of all states.
    JSdiv_aggregated_hmm(i) = total_hmm / 5;      % (reps x conditions) = (3 x 9).
    JSdiv_aggregated_hsmm(i) = total_hsmm / 5;

end

corrected_significance_means = compare_groups(JSdiv_aggregated_hmm, JSdiv_aggregated_hsmm)

average_JSdiv_states_hmm = sum(reshape(JSdiv_states_hmm, 3,3,3,5),3) / 3;
average_JSdiv_states_hsmm = sum(reshape(JSdiv_states_hsmm, 3,3,3,5),3) / 3;
figure
for i = 1:5
   subplot(2,3,i)
   heatmap(average_JSdiv_states_hmm(:,:,:,i));
   caxis([0,0.4])
   colorbar('off')
   title(strcat('State ', num2str(i)))
end
sgtitle('State-specific Jensen-Shannon Divergence, HMM')

figure
for i = 1:5
   subplot(2,3,i)
   heatmap(average_JSdiv_states_hsmm(:,:,:,i));
   caxis([0,0.4])
   colorbar('off')
   title(strcat('State ', num2str(i)))
end
sgtitle('State-specific Jensen-Shannon Divergence, HSMM')

JSdiv_aggregated_hmm = sum(JSdiv_aggregated_hmm,3) / 3;
JSdiv_aggregated_hsmm = sum(JSdiv_aggregated_hsmm,3) / 3;

figure
% Plot aggregated JSdiv
    % Axis
    STN = [2 5 10];
    Coupling = [0.3 0.5 0.8];
    Z_hmm = JSdiv_aggregated_hmm;
    Z_hsmm = JSdiv_aggregated_hsmm;
    
    % Heatmap HMM
    subplot(1,2,1)
    heatmap(STN, Coupling, Z_hmm);
    caxis([0,0.15])
    colorbar('off');
    title('HMM')
    xlabel('STN')
    ylabel('Coupling')
    
    % Heatmap HSMM
    subplot(1,2,2)
    heatmap(STN, Coupling, Z_hsmm);
    caxis([0,0.15])
    colorbar('off');
    title('HSMM')
    xlabel('STN')
    ylabel('Coupling')

sgtitle('Jensen-Shannon divergence (aggregated over states)');

   %{
figure
for i = 1:9
    % Plot individual states JSDiv
    for s = 1:5
        
        subplot(5,2,2*s -1)
        Z_hmm = JSdiv_states_hmm(:,:,s);
        %surf(X,Y,Z_hmm)
        heatmap(Coupling, STN, Z_hmm);
        colorbar('off');
        title(strcat('HMM, state ', num2str(s)))
        xlabel('Coupling')
        ylabel('STN')
        %zlabel('JSdiv')
        
        
        subplot(5,2,2*s)
        Z_hsmm = JSdiv_states_hsmm(:,:,s);
        %surf(X,Y,Z_hsmm)
        heatmap(Coupling, STN, Z_hsmm);
        colorbar('off');
        title(strcat('HSMM, state ', num2str(s)))
        xlabel('Coupling')
        ylabel('STN')
        %zlabel('JSdiv')
        
    end
end
sgtitle('Decoding JSdiv (per state)');
%}
