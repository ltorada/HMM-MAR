function s = nl_logsumexp(a, dim)

    % logsumexp trick for dimension dim of matrix/vector a.
    % Source: https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/39653/versions/1/previews/Shuang_EP_Distributed/utils/logsumexp.m/index.html
    
    if nargin < 2
        dim = 1;
    end
    
    % subtract the largest in each column
    [y, ~] = max(-a,[],dim);
    dims = ones(1,ndims(a));
    dims(dim) = size(a,dim);
    a = -a - repmat(y, dims);
    % Back to neglogs
    s = -y - log(sum(exp(a),dim));
    m = ~isfinite(y);
    if ~isempty(m)
        s(m) = y(m);
    end

end
