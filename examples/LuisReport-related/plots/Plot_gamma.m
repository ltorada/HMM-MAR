v0 = results(i).hsmmtrue.states_seq;
Gammatrue0 = zeros(results(i).timepoints, results(i).hsmmtrue.K);
for t = 1:results(i).timepoints
    Gammatrue0(t,v0(t)) = 1;
end

g = results(i).HMM.gamma;
g2 = results(i).HSMM.gamma;

% Convert Viterbi paths' labels to match the true state labels (match states).
matching_hmm = results(i).HMM.true_to_inferred_mapping;
matching_hsmm = results(i).HSMM.true_to_inferred_mapping;
    
% Get state indexes (excluding repetitions) in zig-zag order left-rigth and top-bottom.
    indexes_HMM = [];
    redundant = [];
    for row = 1:results(i).hsmmtrue.K
       index = find(results(i).HMM.true_to_inferred_mapping(row,:)); % get row.
       new_index = index(find(~ismember(index, indexes_HMM))); % get those not included before from the row.
       if length(new_index) > 1
           redundant = [redundant, new_index(2:end)];
       end
       if ~isempty(new_index)
           indexes_HMM = [indexes_HMM, new_index(1)];
       end
    end
    indexes_HMM = [indexes_HMM, index(find(~ismember(index, indexes_HMM)))]; % include redundant states if not included before.
    
    indexes_HSMM = [];
    for row = 1:results(i).hsmmtrue.K
       index = find(results(i).HSMM.true_to_inferred_mapping(row,:));
       new_index = index(find(~ismember(index, indexes_HSMM)));
       if length(new_index) > 1
           redundant = [redundant, new_index(2:end)];
       end
       if ~isempty(new_index)
           indexes_HSMM = [indexes_HSMM, new_index(1)];
       end
    end
    indexes_HSMM = [indexes_HSMM, index(find(~ismember(index, indexes_HSMM)))];
        
    
figure 

subplot(3,1,1)
bar(Gammatrue0(1:2000,:), 'stacked')
title('Truth')
xlabel('p')
ylabel('Time')
legend('State 1', 'State 2', 'State 3', 'State 4', 'State 5');

subplot(3,1,2)
bar(g(1:2000,indexes_HMM), 'stacked')
title('HMM')
xlabel('p')
ylabel('Time')

subplot(3,1,3)
bar(g2(1:2000,indexes_HSMM), 'stacked')
title('HSMM')
xlabel('p')
ylabel('Time')

sgtitle(strcat(...,
           strcat(...,
              strcat(...,
                 strcat(...,
                 'State Timecourse Posterior. (STN = ', num2str(results(i).STN)),...,
              ', Coupling = '),...
            num2str(results(i).Coupling)),...,
         ')'));