STN = [2 2 2 5 5 5 10 10 10];
Coupling = [0.3 0.5 0.8 0.3 0.5 0.8 0.3 0.5 0.8];

for i = 1:9
    
    figure;
    
    % Get state indexes (excluding repetitions) in zig-zag order left-rigth and top-bottom.
    indexes_HMM = [];
    for row = 1:results(i).hsmmtrue.K
       index = find(results(i).HMM.true_to_inferred_mapping(row,:));
       new_index = index(find(~ismember(index, indexes_HMM)));
       indexes_HMM = [indexes_HMM, new_index];
    end
    
    indexes_HSMM = [];
    for row = 1:results(i).hsmmtrue.K
       index = find(results(i).HSMM.true_to_inferred_mapping(row,:));
       new_index = index(find(~ismember(index, indexes_HSMM)));
       indexes_HSMM = [indexes_HSMM, new_index];
    end
    
    
    % True expected lifetimes
    if results(i).hsmmtrue.sojourns_form == "Geometric"
        means = a/results(i).hsmmtrue.sojourns_parameters;
        subplot(3,2,1)
        for s = 1:hsmmtrue.K
            plot(geopdf(1:results(i).hsmmtrue.M, results(i).hsmmtrue.sojourns_parameters(s)))
            hold on
        end
    else
        means = results(i).hsmmtrue.sojourns_parameters(:,1) .* results(i).hsmmtrue.sojourns_parameters(:,2);
        subplot(3,2,1)
        for s = 1:results(i).hsmmtrue.K
            plot(gampdf(1:results(i).hsmmtrue.M, results(i).hsmmtrue.sojourns_parameters(s,1), results(i).hsmmtrue.sojourns_parameters(s,2)))
            hold on
        end
    end
    
    
    title('Data expected lifetime distributions')
    xlabel('lifetime, d')
    ylabel('p(d)')
    leg = legend(strcat(strcat(strcat('state ',num2str(1), ': mean = '), num2str(means(1)))),...
        strcat(strcat(strcat('state ',num2str(2), ': mean = '), num2str(means(2)))),...
        strcat(strcat(strcat('state ',num2str(3), ': mean = '), num2str(means(3)))),...
        strcat(strcat(strcat('state ',num2str(4), ': mean = '), num2str(means(4)))),...
        strcat(strcat(strcat('state ',num2str(5), ': mean = '), num2str(means(5)))));
    
    
    % True empirical lifetimes
    ss = get_empirical_lifetimes(results(i).hsmmtrue.states_seq, results(i).hsmmtrue.K, results(i).hsmmtrue.M);
    means = round(get_mean_lifetime(ss),1);
    [~, modes] = max(ss);
    
    subplot(3,2,2)
    plot(ss)
    title('Data empirical lifetime distributions')
    xlabel('lifetime, d')
    ylabel('p(d)')
    leg = legend(strcat(strcat(strcat(strcat(strcat('state ',num2str(1), ': mode = '), num2str(modes(1))), ', mean = '), num2str(means(1)))),...
        strcat(strcat(strcat(strcat(strcat('state ',num2str(2), ': mode = '), num2str(modes(2))), ', mean = '), num2str(means(2)))),...
        strcat(strcat(strcat(strcat(strcat('state ',num2str(3), ': mode = '), num2str(modes(3))), ', mean = '), num2str(means(3)))),...
        strcat(strcat(strcat(strcat(strcat('state ',num2str(4), ': mode = '), num2str(modes(4))), ', mean = '), num2str(means(4)))),...
        strcat(strcat(strcat(strcat(strcat('state ',num2str(5), ': mode = '), num2str(modes(5))), ', mean = '), num2str(means(5)))));
    
     % HMM modelled lifetimes
     geo_pmf = zeros(results(i).hsmmtrue.M, results(i).hsmmtrue.K);
     geo_params = diag(results(1).HMM.model.P);
     for s = 1:results(i).hsmmtrue.K
         geo_pmf(:,s) = geopdf(1:100, geo_params(s));
     end
     
    subplot(3,2,3)
    for state = indexes_HMM
        plot(geo_pmf(:,state))
        hold on;
    end
    %[~, modes] = max(results(i).HMM.pd);
    title('Modelled lifetime distributions, HMM')
    xlabel('lifetime, d')
    ylabel('p(d)')
    means = round(get_mean_lifetime(geo_pmf),1);
    [~, modes] = max(geo_pmf);
    leg = legend(strcat(strcat(strcat(strcat(strcat('state ',num2str(indexes_HMM(1)), ': mode = '), num2str(modes(indexes_HMM(1))))), ', mean = '), num2str(means(indexes_HMM(1)))),...
        strcat(strcat(strcat(strcat(strcat('state ',num2str(indexes_HMM(2)), ': mode = '), num2str(modes(indexes_HMM(2))))), ', mean = '), num2str(means(indexes_HMM(2)))),...
        strcat(strcat(strcat(strcat(strcat('state ',num2str(indexes_HMM(3)), ': mode = '), num2str(modes(indexes_HMM(3))))), ', mean = '), num2str(means(indexes_HMM(3)))),...
        strcat(strcat(strcat(strcat(strcat('state ',num2str(indexes_HMM(4)), ': mode = '), num2str(modes(indexes_HMM(4))))), ', mean = '), num2str(means(indexes_HMM(4)))),...
        strcat(strcat(strcat(strcat(strcat('state ',num2str(indexes_HMM(5)), ': mode = '), num2str(modes(indexes_HMM(5))))), ', mean = '), num2str(means(indexes_HMM(5)))));
    
    
     % HMM empirical lifetimes
    subplot(3,2,4)
    for state = indexes_HMM
        plot(results(i).HMM.pd(:,state))
        hold on;
    end
    %[~, modes] = max(results(i).HMM.pd);
    title('Empirical lifetime distributions, HMM')
    xlabel('lifetime, d')
    ylabel('p(d)')
    means = round(get_mean_lifetime(results(i).HMM.pd),1);
    [~, modes] = max(results(i).HMM.pd);
    leg = legend(strcat(strcat(strcat(strcat(strcat('state ',num2str(indexes_HMM(1)), ': mode = '), num2str(modes(indexes_HMM(1))))), ', mean = '), num2str(means(indexes_HMM(1)))),...
        strcat(strcat(strcat(strcat(strcat('state ',num2str(indexes_HMM(2)), ': mode = '), num2str(modes(indexes_HMM(2))))), ', mean = '), num2str(means(indexes_HMM(2)))),...
        strcat(strcat(strcat(strcat(strcat('state ',num2str(indexes_HMM(3)), ': mode = '), num2str(modes(indexes_HMM(3))))), ', mean = '), num2str(means(indexes_HMM(3)))),...
        strcat(strcat(strcat(strcat(strcat('state ',num2str(indexes_HMM(4)), ': mode = '), num2str(modes(indexes_HMM(4))))), ', mean = '), num2str(means(indexes_HMM(4)))),...
        strcat(strcat(strcat(strcat(strcat('state ',num2str(indexes_HMM(5)), ': mode = '), num2str(modes(indexes_HMM(5))))), ', mean = '), num2str(means(indexes_HMM(5)))));
    
    
    % HSMM expected lifetimes
    subplot(3,2,5)
    for state = indexes_HSMM
        plot(results(i).HSMM.model.sojourns(:,state))
        hold on;
    end
    title('Modelled lifetime distributions, HSMM')
    xlabel('lifetime, d')
    ylabel('p(d)')
    %means = round(get_mean_lifetime(results(i).HSMM.model.sojourns),1);
    means = round(results(i).HSMM.model.sojourn_means ,1);
    %[~, modes] = max(results(i).HSMM.model.sojourns);
    modes = results(i).HSMM.model.sojourn_modes_before_smoothing;
    leg = legend(strcat(strcat(strcat(strcat(strcat('state ',num2str(indexes_HSMM(1)), ': mode = '), num2str(modes(indexes_HSMM(1))))), ', mean = '), num2str(means(indexes_HSMM(1)))),...
        strcat(strcat(strcat(strcat(strcat('state ',num2str(indexes_HSMM(2)), ': mode = '), num2str(modes(indexes_HSMM(2))))), ', mean = '), num2str(means(indexes_HSMM(2)))),...
        strcat(strcat(strcat(strcat(strcat('state ',num2str(indexes_HSMM(3)), ': mode = '), num2str(modes(indexes_HSMM(3))))), ', mean = '), num2str(means(indexes_HSMM(3)))),...
        strcat(strcat(strcat(strcat(strcat('state ',num2str(indexes_HSMM(4)), ': mode = '), num2str(modes(indexes_HSMM(4))))), ', mean = '), num2str(means(indexes_HSMM(4)))),...
        strcat(strcat(strcat(strcat(strcat('state ',num2str(indexes_HSMM(5)), ': mode = '), num2str(modes(indexes_HSMM(5))))), ', mean = '), num2str(means(indexes_HSMM(5)))));
    
   
    % HSMM empirical lifetimes
    subplot(3,2,6)
    for state = indexes_HSMM
        plot(results(i).HSMM.pd(:,state))
        hold on;
    end
    title('Empirical lifetime distributions, HSMM')
    xlabel('lifetime, d')
    ylabel('p(d)')
    [~, modes] = max(results(i).HSMM.pd);
    means = round(get_mean_lifetime(results(i).HSMM.pd),1);
    leg = legend(strcat(strcat(strcat(strcat(strcat('state ',num2str(indexes_HSMM(1)), ': mode = '), num2str(modes(indexes_HSMM(1))))), ', mean = '), num2str(means(indexes_HSMM(1)))),...
        strcat(strcat(strcat(strcat(strcat('state ',num2str(indexes_HSMM(2)), ': mode = '), num2str(modes(indexes_HSMM(2))))), ', mean = '), num2str(means(indexes_HSMM(2)))),...
        strcat(strcat(strcat(strcat(strcat('state ',num2str(indexes_HSMM(3)), ': mode = '), num2str(modes(indexes_HSMM(3))))), ', mean = '), num2str(means(indexes_HSMM(3)))),...
        strcat(strcat(strcat(strcat(strcat('state ',num2str(indexes_HSMM(4)), ': mode = '), num2str(modes(indexes_HSMM(4))))), ', mean = '), num2str(means(indexes_HSMM(4)))),...
        strcat(strcat(strcat(strcat(strcat('state ',num2str(indexes_HSMM(5)), ': mode = '), num2str(modes(indexes_HSMM(5))))), ', mean = '), num2str(means(indexes_HSMM(5)))));
    
    sgtitle(strcat(strcat(strcat(strcat('STN = ', num2str(STN(i))), ', Coupling = '), num2str(Coupling(i))), strcat(';  timepoints = ', num2str(results(i).timepoints))))
    
end