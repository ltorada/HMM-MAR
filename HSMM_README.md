## Training a HSMM

### 3 new options
options.M:                             Maximal lifetime duration allowed (default to 100 timepoints) [int].\
options.hsmm:                          Whether you want to train a hsmm or not (in which case a hmm will be trained) [1/0].\
options.sojourns_init (optional):      Lifetime distribution with which to initialise the hsmm [(M x S) matrix].\

## Example

```
addpath(genpath('.'))

STN = 5;
Coupling = 0.8; 
timepoints = 10000;\
sim_model_definition;\                                              # Create hsmmtrue structure (used for the simulation, or simply by run_hmmar if at contains at least the fields M, T and K).

[X_hsmm, Gammatrue_hsmm, states_seq] = simhsmm(hsmmtrue);           # Simulate data and keep true state timecourse.

figure
subplot(2,1,1)
plot(X_hsmm(1:100,:))                                               # Plot all simualted dimensions during 100 timepoints.
subplot(2,1,2)
bar(Gammatrue_hsmm(1:100, :), 'stacked')                            # Plot states timecourse during 100 timepoints.
title('HSMM Viterbi path')
ylabel('p')
xlabel('Time')

[hsmm, gamma, ~, vpath] = run_hmmar(X_hsmm, hsmmtrue, "hsmm", 0);   # Train hsmm.

figure
bar(gamma(1:100, :), 'stacked')                                      # Plot states timecourse during 100 timepoints.
title('HSMM states posterior')
ylabel('p')
xlabel('Time')
```

## Modifications with respect to the only-HMM version of OHBA HMM-MAR:

### hmmtrain.m
hsmm version of the Viterbi used if options.hsmm was set to 1.


### hsinference.m
hsmm version of the Forward-Backwards used if options.hsmm was set to 1.


### hmmmar_init.m (stochastic version not adapted yet)
Initialisation of the state timecourse with a semi-Markov chain (with uniform lifetime distributions between 1 and M) if options.hsmm was set to 1.

### Still needs to be adapted

'eval' folder (to compute the free-energy between q(H) and p(H), which now differs in the unnormalized posterior p(S,X)).
