Accuracy_aggregated_hmm = zeros(3,3,3); % (Coupling x STN)
Accuracy_aggregated_hsmm = zeros(3,3,3); % (Coupling x STN)

Accuracy_states_hmm = zeros(27,5);
Accuracy_states_hsmm = zeros(27,5);

%ref_switches = zeros(1,27);
hmm_redundant = zeros(1,27);
hsmm_redundant = zeros(1,27);

for i = 1:27
    
    
    matching_hmm = results(i).HMM.true_to_inferred_mapping;
    matching_hsmm = results(i).HSMM.true_to_inferred_mapping;
    
    hmm_redundant(i) = sum(sum(matching_hmm,2) > 1);
    hsmm_redundant(i) = sum(sum(matching_hsmm,2) > 1);

    seq_hmm = results(i).HMM.vpath';
    seq_hsmm = results(i).HSMM.vpath;
    reference = results(i).hsmmtrue.states_seq;
    
    %hmm_matches = zeros(1,5);
    %hsmm_matches = zeros(1,5);
    
    % Get accuracy of each state under this condition (i).
    for s = 1:5
        % HMM  
            candidates = [];
            hmm_ms = find(matching_hmm(s,:));
            for ms = hmm_ms % choose best accuracy only, assuming rendundancy in the rest and therefore penalizing for it.
                candidates = [candidates, sum(((seq_hmm == ms) + (reference == s)) == 2)];
            end
            [Accuracy_states_hmm(i,s), ~] = max(candidates / sum(reference == s));
            
            %for ms = find(matching_hmm(s,:))
            %    hmm_matches(s) = hmm_matches(s) + sum(((seq_hmm == ms) + (reference == s)) == 2);
            %end
            %Accuracy_states_hmm(i,s) = hmm_matches(s) / sum(reference == s);
            
        % HSMM
        
            candidates = [];
            hsmm_ms = find(matching_hsmm(s,:));
            for ms = hsmm_ms % choose best accuracy only, assuming rendundancy in the rest and therefore penalizing for it.
                candidates = [candidates, sum(((seq_hsmm == ms) + (reference == s)) == 2)];
            end
            [Accuracy_states_hsmm(i,s), ~] = max(candidates / sum(reference == s));
            
            %for ms = find(matching_hsmm(s,:))
            %    hsmm_matches(s) = hsmm_matches(s) + sum(((seq_hsmm == ms) + (reference == s)) == 2);
            %end
            %Accuracy_states_hsmm(i,s) = hsmm_matches(s) / sum(reference == s);
    end
    
    
    %Accuracy_aggregated_hmm(i) = sum(hmm_matches) / length(reference);
    %Accuracy_aggregated_hsmm(i) = sum(hsmm_matches) / length(reference);
    
    %ref_switches(i) = sum(diff(reference)~=0);
    %hmm_switches(i) = sum(diff(seq_hmm)~=0);
    %hsmm_switches(i) = sum(diff(seq_hsmm)~=0);
end

hmm_redundant = sum(reshape(hmm_redundant, 3,3,3),3)
hsmm_redundant = sum(reshape(hsmm_redundant, 3,3,3),3)

Average_accuracy_per_state_hmm = sum(reshape(Accuracy_states_hmm, 3,3,3,5),3) / 3;
Average_accuracy_per_state_hsmm = sum(reshape(Accuracy_states_hsmm, 3,3,3,5),3) / 3;
figure
for i = 1:5
   subplot(2,3,i)
   heatmap(Average_accuracy_per_state_hmm(:,:,:,i));
   colorbar('off')
   title(strcat('State ', num2str(i)))
end
sgtitle('State-specific Decoding Accuracy HMM')

figure
for i = 1:5
   subplot(2,3,i)
   heatmap(Average_accuracy_per_state_hsmm(:,:,:,i));
   colorbar('off')
   title(strcat('State ', num2str(i)))
end
sgtitle('State-specific Decoding Accuracy HSMM')

Accuracy_aggregated_hmm = sum(reshape(sum(Accuracy_states_hmm,2), 3,3,3), 3) / 15;
Accuracy_aggregated_hsmm = sum(reshape(sum(Accuracy_states_hsmm,2), 3,3,3), 3) / 15;

%Accuracy_aggregated_hmm = sum(Accuracy_aggregated_hmm,3) / 3;
%Accuracy_aggregated_hsmm = sum(Accuracy_aggregated_hsmm,3) / 3;
    
    %{
    % Average-out all the states under this condition (i).
    avg = sum(Accuracy_states_hmm, 3) / 5;
    Accuracy_aggregated_hmm(i) = avg(i);
    avg = sum(Accuracy_states_hsmm, 3) / 5;
    Accuracy_aggregated_hsmm(i) = avg(i);

    % Get number of switches per state.
    ref_switches = zeros(1,5);
    hmm_switches = zeros(1,5);
    hsmm_switches = zeros(1,5);
    for s = 1:5
        ref_switches(s) = sum(diff(reference == s)~=0);
        hmm_switches(s) = sum(diff(m_seq_hmm == s)~=0);
        hsmm_switches(s) = sum(diff(m_seq_hsmm == s)~=0);
    end
    aggregated_ref_switches = sum(ref_switches);
    aggregated_hmm_switches = sum(hmm_switches);
    aggregated_hsmm_switches = sum(hsmm_switches);
    
    display_aggregated_matrix = [ref_switches; hmm_switches; hsmm_switches];
    display_states_matrix = [aggregated_ref_switches; aggregated_hmm_switches; aggregated_hsmm_switches];
    %}
    % Plot aggregated accuracy
    STN = [2 5 10];
    Coupling = [0.3 0.5 0.8];
    Z_hmm = Accuracy_aggregated_hmm;
    Z_hsmm = Accuracy_aggregated_hsmm;
    
    figure
    
    % Plot avera accuracy for the HMM.
    subplot(1,2,1)
    heatmap(STN,Coupling, Z_hmm);
    caxis([0,1])
    colorbar('off')
    %title(strcat('HMM. Number of switches:  ', num2str(hmm_switches)))
    title('HMM')
    ylabel('Coupling')
    xlabel('STN')
    
    % Plot avera accuracy for the HSMM.
    subplot(1,2,2)
    heatmap(STN, Coupling, Z_hsmm);
    caxis([0,1])
    colorbar('off')
    %title(strcat('HSMM. Number of switches:   ', num2str(hsmm_switches)))
    title('HSMM')
    ylabel('Coupling')
    xlabel('STN')
    
sgtitle('Decoding Accuracy (aggregated over states)');

%{
figure
diff_hmm = hmm_switches - ref_switches;
diff_hsmm = hsmm_switches - ref_switches;
subplot(1,2,1)
heatmap(diff_hmm)
subplot(1,2,2)
heatmap(diff_hsmm)
%}
    %{
    figure
for i = 1:9
    % Plot individual states accuracy
    for s = 1:5
        
        subplot(5,2,2*s -1)
        Z_hmm = Accuracy_states_hmm(:,:,s);
        %surf(X,Y,Z_hmm)
        heatmap(Coupling, STN, Z_hmm);
        colorbar('off')
        title(strcat('HMM, state ', num2str(s)))
        xlabel('Coupling')
        ylabel('STN')
        %zlabel('DA')
        
        
        subplot(5,2,2*s)
        Z_hsmm = Accuracy_states_hsmm(:,:,s);
        %surf(X,Y,Z_hsmm)
        heatmap(Coupling, STN, Z_hsmm);
        colorbar('off')
        title(strcat('HSMM, state ', num2str(s)))
        xlabel('Coupling')
        ylabel('STN')
        %zlabel('DA')
        
    end

end
    sgtitle('Decoding Accuracy (per state)');
%}