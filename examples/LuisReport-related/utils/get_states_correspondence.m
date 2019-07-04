function true_to_inferred_mapping = get_states_correspondence(results, hsmmtrue)

% Get Boolean matrix mapping true to inferred states.
%
% INPUTS
% ------
%
% results     structure corresponding to the trained model.
% hsmmtrue    structure corresponding to the true model (used for simulating the data).
%
%
% OUTPUT
% ------
% 
% true_to_inferred_mapping     Boolean matrix
%
%
% Author: Luis Torada Aguilella
    
    
    states_correspondence = zeros(hsmmtrue.K, 2);
    reverse_states_correspondence = zeros(hsmmtrue.K, 2);
    
    % true state --> inferred state
    difs = zeros(1,hsmmtrue.K);
    for s = 1:hsmmtrue.K
        true_cov = hsmmtrue.state(s).Cov;                        % Get covariance matrix for true state s.
        for s2 = 1:hsmmtrue.K                                    % Find most simular one from the trained model.
            [inferred_cov, ~] = getFuncConn(results.model, s2);
            % Normalize with diagonal elements (variances).
            norm_inferred_cov = tril(inferred_cov ./ diag(inferred_cov));
            norm_inferred_cov = norm_inferred_cov + norm_inferred_cov' - eye(hsmmtrue.K);
            norm_true_cov = tril(true_cov ./ diag(true_cov));
            norm_true_cov = norm_true_cov + norm_true_cov' - eye(hsmmtrue.K);
            % Calculate sum(residuals) between them.
            difs(s2) = sum(abs(norm_inferred_cov - norm_true_cov), 'all');
         end
        [~, ind] = min(difs); % State from HMM corresponding to simualted state s.
        states_correspondence(s, :) = [s, ind];
    end
    
    % inferred state --> true state   (same as before but in the reverse order)
    difs = zeros(1,hsmmtrue.K);
    for s2 = 1:hsmmtrue.K
        [inferred_cov, ~] = getFuncConn(results.model, s2);
        for s = 1:hsmmtrue.K
            true_cov = hsmmtrue.state(s).Cov;
            % Normalize with diagonal elements (variances).
            norm_inferred_cov = tril(inferred_cov ./ diag(inferred_cov));
            norm_inferred_cov = norm_inferred_cov + norm_inferred_cov' - eye(hsmmtrue.K);
            norm_true_cov = tril(true_cov ./ diag(true_cov));
            norm_true_cov = norm_true_cov + norm_true_cov' - eye(hsmmtrue.K);
            % Calculate sum(residuals) between them.
            difs(s) = sum(abs(norm_inferred_cov - norm_true_cov), 'all');
        end
        [~, ind] = min(difs); % State from HMM corresponding to simualted state s.
        reverse_states_correspondence(s2, :) = [s2, ind];
    end
    
    %%% Condense into a matrix (S x S) mapping true states (row indexes) to
    %%% inferred states modelling them (one-hot encoding in the corresponding row).
    true_to_inferred_mapping = zeros(hsmmtrue.K);
    % First, add best approximation to each true state
    for true = 1:hsmmtrue.K
        best_approx = states_correspondence(true,2);
        true_to_inferred_mapping(true,best_approx) = 1;
    end
    % Then, add extra states modelling the true states.
    for inferred = 1:hsmmtrue.K
        if ~ismember(inferred, states_correspondence(:,2))
            best_approx = reverse_states_correspondence(inferred,2);
            true_to_inferred_mapping(best_approx, inferred) = 1;
        end
    end
    
end