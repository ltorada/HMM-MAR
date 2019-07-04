mean_residuals_hmm = zeros(3,3,3);
mean_residuals_hsmm = zeros(3,3,3);

for i = 1:27
%{
    % Get state indexes (excluding repetitions) in zig-zag order left-rigth and top-bottom.
    indexes_HMM = [];
    for row = 1:hsmmtrue.K
       index = find(results(i).HMM.true_to_inferred_mapping(row,:));
       new_index = index(find(~ismember(index, indexes_HMM)));
       indexes_HMM = [indexes_HMM, new_index];
    end
    
    indexes_HSMM = [];
    for row = 1:hsmmtrue.K
       index = find(results(i).HSMM.true_to_inferred_mapping(row,:));
       new_index = index(find(~ismember(index, indexes_HSMM)));
       indexes_HSMM = [indexes_HSMM, new_index];
    end
    %}
    %{
    figure; % One figure per condition, showing one heatmap per model (cols) and state (rows)
    
    % True Cov
    counter = 1;
    for s = 1:hsmmtrue.K
        subplot(3,5,s);
        covmat = hsmmtrue.state(s).Cov;
        heatmap(covmat);
        %colorbar('off');
        title(strcat('Truth. State ', num2str(counter)));
        counter = counter + 1;
    end
    
    % HMM Cov
    counter = 1;
    for s = indexes_HMM
        subplot(3,5,counter + 5);
        [covmat,corrmat] = getFuncConn(results(i).HMM.model, s);
        heatmap(covmat);
        %colorbar('off');
        title(strcat('HMM. State ', num2str(s)));
        counter = counter + 1;
    end
    %}
    
    if i < 10
        rep = 1;
    elseif i < 19
        rep = 2;
    else
        rep = 3;
    end
    
    % Get residuals
    total_residuals = 0;
    for s = 1:5
        mss = find(results(i).HMM.true_to_inferred_mapping(s,:));
        total_residuals_state = 0;
        for ms = mss
            [covmat,~] = getFuncConn(results(i).HMM.model, ms);
            total_residuals_state = total_residuals_state + sum(abs(covmat - hsmmtrue.state(s).Cov), 'all'); 
        end
        total_residuals_state = total_residuals_state / length(mss);
        total_residuals = total_residuals + total_residuals_state; % Add total residuals of each state.
    end
    mean_residuals_hmm(i) = total_residuals / (5*5*5); % In each condition, take the mean over all states (5) and cov matrix dimensions (5*5) == 5*5*5;

    %{
    % HSMM Cov
    counter = 1;
    for s = indexes_HSMM
        subplot(3,5,counter + 10);
        [covmat,corrmat] = getFuncConn(results(i).HSMM.model, s);
        heatmap(covmat);
        %colorbar('off');
        title(strcat('HSMM. State ', num2str(s)));
        counter = counter + 1;
    end
    %}
    total_residuals = 0;
    for s = 1:5
        mss = find(results(i).HSMM.true_to_inferred_mapping(s,:));
        total_residuals_state = 0;
        for ms = mss
            [covmat,~] = getFuncConn(results(i).HSMM.model, ms);
            total_residuals_state = total_residuals_state + sum(abs(covmat - hsmmtrue.state(s).Cov), 'all');
        end
        total_residuals_state = total_residuals_state / length(mss);
        total_residuals = total_residuals + total_residuals_state;
    end
    mean_residuals_hsmm(i) = total_residuals / (5*5*5);

    %sgtitle(strcat( strcat(strcat(strcat(strcat('STN = ', num2str(results(i).STN)), ', Coupling = '), num2str(results(i).Coupling)), ';  timepoints = '), num2str(timepoints)));
    
end

mean_residuals_hmm = sum(mean_residuals_hmm,3);
mean_residuals_hsmm = sum(mean_residuals_hsmm,3);

%mean_residuals_hmm = reshape(mean_residuals_hmm, 3,3);
%mean_residuals_hsmm = reshape(mean_residuals_hsmm, 3,3);

STN = [2 5 10];
Coupling = [0.3 0.5 0.8];
    
figure
subplot(1,2,1)
heatmap(STN, Coupling, mean_residuals_hmm);
caxis([0 0.65])
colorbar('off')
xlabel('STN')
ylabel('Coupling')
title('HMM')

subplot(1,2,2)
heatmap(STN, Coupling, mean_residuals_hsmm);
caxis([0 0.65])
colorbar('off')
xlabel('STN')
ylabel('Coupling')
title('HSMM')

sgtitle('Mean covariance residuals')

