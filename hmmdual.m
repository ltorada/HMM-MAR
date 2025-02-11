function [hmm,Gamma] = hmmdual(data,T,hmm,Gamma,Xi,residuals)
%
% Dual estimation of the HMM
%
% INPUTS:
%
% data          observations - a struct with X (time series) and C (classes; optional)
% T             Number of time points for each time series
% hmm           hmm structure with options specified in hmm.train
% Gamma         Initial state courses
% Xi            joint probability of past and future states conditioned on data
% residuals     in case we train on residuals, the value of those.
%
% OUTPUTS
% hmm           estimated HMMMAR model
% Gamma         estimated p(state | data)
%
% Author: Diego Vidaurre, OHBA, University of Oxford (2019)

% to fix potential compatibility issues with previous versions
hmm = versCompatibilityFix(hmm);

if nargin<6, residuals = []; end
if nargin<5, Xi = []; end
if nargin<4, Gamma = []; end

if iscell(T)
    for i = 1:length(T)
        if size(T{i},1)==1, T{i} = T{i}'; end
    end
    if size(T,1)==1, T = T'; end
    T = cell2mat(T);
end
checkdatacell;
N = length(T);

train = hmm.train;
checkdatacell;
data = data2struct(data,T,train);
% Filtering
if ~isempty(train.filter)
    data = filterdata(data,T,train.Fs,train.filter);
end
% Detrend data
if train.detrend
    data = detrenddata(data,T);
end
% Standardise data and control for ackward trials
data = standardisedata(data,T,train.standardise);
% Leakage correction
if train.leakagecorr ~= 0
    data = leakcorr(data,T,train.leakagecorr);
end
% Hilbert envelope
if train.onpower
    data = rawsignal2power(data,T);
end
% Leading Phase Eigenvectors
if train.leida
    data = leadingPhEigenvector(data,T);
end
% pre-embedded  PCA transform
if length(train.pca_spatial) > 1 || train.pca_spatial > 0
    if isfield(train,'As')
        data.X = bsxfun(@minus,data.X,mean(data.X));
        data.X = data.X * train.As;
    else
        [train.As,data.X] = highdim_pca(data.X,T,train.pca_spatial);
    end
end
% Embedding
if length(train.embeddedlags) > 1
    [data,T] = embeddata(data,T,train.embeddedlags);
end
% PCA transform
if length(train.pca) > 1 || train.pca > 0
    if isfield(train,'A')
        data.X = bsxfun(@minus,data.X,mean(data.X));
        data.X = data.X * train.A;
    else
        [train.A,data.X] = highdim_pca(data.X,T,train.pca,0,0,0,train.varimax);
    end
    % Standardise principal components and control for ackward trials
    data = standardisedata(data,T,train.standardise_pc);
    train.ndim = size(train.A,2);
    train.S = ones(train.ndim);
    orders = formorders(train.order,train.orderoffset,train.timelag,train.exptimelag);
    train.Sind = formindexes(orders,train.S);
end
% Downsampling
if train.downsample > 0
    [data,T] = downsampledata(data,T,train.downsample,train.Fs);
end

if isempty(residuals)
    if ~isfield(hmm.train,'Sind')
        orders = formorders(hmm.train.order,hmm.train.orderoffset,hmm.train.timelag,hmm.train.exptimelag);
        hmm.train.Sind = formindexes(orders,hmm.train.S);
    end
    residuals =  getresiduals(data.X,T,hmm.train.Sind,hmm.train.maxorder,hmm.train.order,...
        hmm.train.orderoffset,hmm.train.timelag,hmm.train.exptimelag,hmm.train.zeromean);
end


if isempty(Gamma) || isempty(Xi)
    [Gamma,~,Xi] = hsinference(data,T,hmm,residuals); 
end
setxx;

hmm = obsupdate(T,Gamma,hmm,residuals,XX,XXGXX);
hmm = hsupdate(Xi,Gamma,T,hmm);

Gamma = hsinference(data,T,hmm,residuals); 

end

