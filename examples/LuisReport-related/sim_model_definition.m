%
% Define the HSMM model that will be used for simulating a dataset.
%
% OUTPUT
% ------
%
% hsmmtrue structure containing all needed for simulating a dataset.
%
%
% Author: Luis Torada Aguilella


% Observations parameters
hsmmtrue = struct();
ndim = 5;                                             % number of channels.
hsmmtrue.T = timepoints;                                              % number of data points.

% Markov chain parameters
hsmmtrue.K = 5;                                       % number of states.
hsmmtrue.P = 1/2 * (ones(hsmmtrue.K) - diag(ones(1,hsmmtrue.K)));       % P MUST BE ZERO DIAGONAL.
hsmmtrue.Pi = 1/3 * ones(1,hsmmtrue.K);

%sojourn_form = "Geometric"; 
%sojourn_parameters = [0.1; 0.2; 0.3; 0.4; 0.5];
%sojourn_parameters = 0.3 * ones(hsmmtrue.K,1);

sojourn_form = "Gamma";                           % Parametric form of the sojourn distributions.                                                      % Alternatives: 
sojourn_parameters = cat(2, [5 20 40 60 100]', 0.5*ones(1,hsmmtrue.K)');
sojourn_parameters = [5 0.5; 10 0.5];

hsmmtrue.M = 100;

noise = 0.1;
if exist("STN")
    activation_noise = STN * noise;
else
    activation_noise = 0.5;
end
if exist("Coupling")
    coupling_noise = Coupling * (activation_noise - noise);
else
    coupling_noise = 0.75 * (activation_noise - noise);
end
sim_mean_differences = 1;

networks = [1 3; 2 4; 3 5; 4 1; 5 2];


% -- Observation model initialization. (Modify this block for advanced control only) -- %
for k = 1:hsmmtrue.K
    hsmmtrue.state(k).W.Mu_W = 0 * sim_mean_differences * ones(1,ndim);
  
    hsmmtrue.state(k).Cov = noise * eye(ndim);
    hsmmtrue.state(k).Cov(networks(k,1), networks(k,1)) = activation_noise;
    hsmmtrue.state(k).Cov(networks(k,2), networks(k,2)) = activation_noise;
    hsmmtrue.state(k).Cov(networks(k,1), networks(k,2)) = coupling_noise;
    hsmmtrue.state(k).Cov( networks(k,2), networks(k,1) ) = coupling_noise;
    
end

% -- Sojourn pmfs initialization. (Modify this block for advanced control only) -- %
   
d = inf(hsmmtrue.M, hsmmtrue.K);

if sojourn_form == "Geometric"
    for j = 1:hsmmtrue.K
           d(:,j) = geopdf(1:hsmmtrue.M, sojourn_parameters(j,:)) ./ sum(geopdf(1:hsmmtrue.M, sojourn_parameters(j,:)));
    end
end
   
if sojourn_form == "Gamma"
    for j = 1:hsmmtrue.K
        d(:,j) = gampdf(1:hsmmtrue.M, sojourn_parameters(j,1), sojourn_parameters(j,2)) ./ sum( gampdf(1:hsmmtrue.M, sojourn_parameters(j,1), sojourn_parameters(j,2)) );
    end
end

%hsmmtrue.sojourns = d;
hsmmtrue.sojourns_form = sojourn_form;
hsmmtrue.sojourns_parameters = sojourn_parameters;

hsmmtrue.train = [];


