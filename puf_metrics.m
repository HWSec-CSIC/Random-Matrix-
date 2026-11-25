function metrics = puf_metrics(C, R, pairing)
%PUF_METRICS Four scalar diagnostics for paired challenge and response data.
%   T = puf_metrics(C,R) evaluates dCor on corresponding entries C(i,j)
%   and R(i,j), as in the grid experiment. C and R must have the same size.
%   T = puf_metrics(C,R,'rows') pairs C(i,:) with R(i,:) as two vectors.
%   Their numbers of features may differ. Rows must be in the same order.
%
%   In both modes R is N observations by P response features. Spectral
%   metrics use all rows and columns of the column-standardized response.
%   rho_eff is the entropy effective rank of the covariance eigenvalues
%   divided by min(N-1,P). S_spike is the summed excess above the fitted
%   Marchenko-Pastur upper edge divided by the total spectral energy.
%   beta_Brody uses Gaussian unfolding and maximum likelihood on [0,1].
%   dCor uses the biased V-statistic with Euclidean distances on the raw
%   data, using all pairs up to 1000 and a fixed common subsample above it.
%   Complex distances include both real and imaginary parts.
%
%   A constant response spectrum of zero energy gives NaN spectral metrics.
%   Brody returns NaN for fewer than 30 usable spacings or degenerate gaps.
%   dCor returns NaN for fewer than three pairs or zero distance variance.

if nargin < 3
    pairing = 'entries';
end

validateattributes(C, {'numeric','logical'}, {'2d','nonempty','finite'});
validateattributes(R, {'numeric','logical'}, {'2d','nonempty','finite'});

[N,P] = size(R);

if N < 2 || size(C,1) ~= N
    error('C and R must have the same number of rows (at least two).');
end

pairing = validatestring(pairing, {'entries','rows'});

if strcmp(pairing,'entries') && ~isequal(size(C),size(R))
    error('Entry pairing requires equal matrix sizes. Use rows for vector CRPs.');
end

C = full(double(C));
R = full(double(R));

% Distance correlation on the original challenge-response pairs.
seed = mod(2718 + 1000003*N + 1000033*P, 2^32-1);

if strcmp(pairing,'entries')
    dCor = distanceCorrelation(C(:),R(:),seed);
else
    dCor = distanceCorrelation(C,R,seed);
end

% Center and standardize each response feature. Keep constant columns.
Z = R - mean(R,1);
Z(:,all(R == R(1,:),1)) = 0;

sd = std(Z,0,1);
active = sd > 0;
Z(:,active) = Z(:,active) ./ sd(active);

% Covariance eigenvalues from the singular values of Z.
sigma = svd(Z,'econ');
lambda = sigma.^2 / (N-1);
total = sum(lambda);

rho_eff = NaN;
S_spike = NaN;
beta_Brody = NaN;

if total > 0
    % Normalized effective rank using spectral energy weights.
    weights = lambda(lambda > 0) / total;
    rho_eff = exp(-sum(weights .* log(weights))) / min(N-1,P);

    % Normalized excess above the upper Marchenko-Pastur edge.
    upperEdge = (total/P) * (1 + sqrt(P/(N-1)))^2;
    S_spike = sum(max(lambda-upperEdge,0)) / total;

    % Brody fit on the numerically nonzero covariance eigenvalues.
    beta_Brody = fitBrody(sort(lambda(sigma > 1e-10*sigma(1))));
end

metrics = table(rho_eff,S_spike,beta_Brody,dCor);

end


function beta = fitBrody(levels)

beta = NaN;
count = numel(levels);

% Discard 10 percent of the levels at each edge before taking spacings.
trim = floor(0.10*count);
interior = (1+trim):(count-trim);

if numel(interior)-1 < 30
    return
end

% Estimate the smooth cumulative level count using all retained levels.
bandwidth = 1.06*std(levels)*count^(-1/5);

if bandwidth <= 0 || ~isfinite(bandwidth)
    return
end

unfolded = zeros(numel(interior),1);

for k = 1:numel(interior)
    offsets = (levels-levels(interior(k))) / (sqrt(2)*bandwidth);
    unfolded(k) = sum(0.5*erfc(offsets));
end

s = diff(unfolded);

if any(s <= 1e-12) || any(~isfinite(s))
    return
end

s = s/mean(s);

% Minimize the negative log-likelihood, including both endpoints.
coef = @(b) exp((b+1)*gammaln((b+2)/(b+1)));
cost = @(b) -sum(log(b+1) + log(coef(b)) + b*log(s) - coef(b)*s.^(b+1));

candidates = [0, fminbnd(cost,0,1), 1];
[~,best] = min(arrayfun(cost,candidates));
beta = candidates(best);

end


function value = distanceCorrelation(X,Y,seed)

value = NaN;
n = size(X,1);

if n < 3
    return
end

% Use the same observations in X and Y, preserving their pairing.
if n > 1000
    stream = RandStream('mt19937ar','Seed',seed);
    idx = sort(randperm(stream,n,1000));
    X = X(idx,:);
    Y = Y(idx,:);
end

A = centeredDistances(X);
B = centeredDistances(Y);

denom = norm(A,'fro') * norm(B,'fro');

if denom > 0
    value = sqrt(max(0,min(1,sum(A(:).*B(:))/denom)));
end

end


function A = centeredDistances(X)

% A common scale factor cancels in the normalized distance correlation.
scale = max(abs(X(:)));

if scale > 0
    X = X / scale;
end

n = size(X,1);

if size(X,2) == 1
    D = abs(X-X.');
else
    D = zeros(n);

    for k = 1:n-1
        delta = X(k+1:end,:) - X(k,:);
        distances = sqrt(sum(abs(delta).^2,2));
        D(k+1:end,k) = distances;
        D(k,k+1:end) = distances.';
    end
end

% Double centering for the biased distance-correlation estimator.
A = D - mean(D,2) - mean(D,1) + mean(D(:));

end