clear; close;
% number of diracs
nDiracs = 2;
% number of moments
nMoments = 5;
% max degree of polynomials
% nDegMax = ceil(2 * nDiracs - 1);
nDegMax = nMoments - 1;
% kernels of finite support
lSignal = 2048;
% sampling period
tSample = 64;
% max amplitue
weightMax = 32;
% number of shifts
nShifts = 31;
% number of iterations
iter = 6;
% time of sampling points
t = 0: 1 / tSample : (lSignal - 1) / tSample;
% variance
sigma = [1e-1; 1e-3; 1e-5];
% noise
noise = sigma * randn(1, nDegMax + 1);
poly = zeros(nDegMax + 1, lSignal);
coeffs = zeros(nDegMax + 1, nShifts + 1);
phi = zeros(1, lSignal);
phiSet = zeros(1 + nShifts, length(phi));
%% Daubechies
% polynomials of max degree N can be reproduced by a scaling function that
% generates wavelets with (N + 1) vanishing moments
[phiT, psiT, xval] = wavefun('dB5', iter);
% extend the scaling function to the signal length
phi(1: length(phiT)) = phiT;
for iShift = 0: nShifts
    % obtain all shifted phi (circular shift)
    phiSet(iShift + 1, :) = circshift(phi, iShift * tSample);
    % obtain all shifted phi (remove tails)
%     phiSet(iShift + 1, :) = [zeros(1, iShift * tSample), phi(1: end - iShift * tSample)];
end
for iDeg = 0: nDegMax
    % all possible polynomials determined by time, max degree nDegMax
    poly(iDeg + 1, :) = t .^ iDeg;
    % coefficients of corresponding kernels
    coeffs(iDeg + 1, :) = dot(repmat(poly(iDeg + 1, :), nShifts + 1, 1), phiSet, 2).' ./ tSample;
end
%% Diracs stream generation
signal = zeros(1, lSignal);
tau = zeros(1, nDegMax + 1);
locations = sort(randperm(lSignal, nDiracs)) / tSample
weights = randperm(weightMax, nDiracs);
signal(locations * tSample) = weights;
samples = signal * phiSet';
for iDeg = 0: nDegMax
    tau(1, iDeg + 1) = dot(coeffs(iDeg + 1, :), samples);
end
tauNoisy = tau + noise(3, :);
tauNoisyMatrix = zeros(nMoments - nDiracs, nDiracs + 1);
for iDirac = 1: nMoments - nDiracs
    tauNoisyMatrix(iDirac, :) = flip(tauNoisy(iDirac: iDirac + nDiracs));
end
[u, lambda, v] = svd(tauNoisyMatrix);
h = v(:, end);
filterCoeffs = h ./ h(1);
flag = 1;
% roots of z-transform of the filter corresponds to the pulse locations:
% H(z) = (1-t0z^(-1))(1-t1z^(-1)) = 1-(t0+t1)z^(-1)+t0t1z(-2)
% h(0)=1; h(1)=t0+t1; h(2)=t0t1
func = @(loc) ([loc(1) * loc(2) - filterCoeffs(3); -(loc(1) + loc(2)) - filterCoeffs(2)]);
% solve to get pulse locations
locationsRec = sort(fsolve(func, [1 1]))
% Vandermonde system
locMatrix = fliplr(vander(locationsRec))';
% weighted sum of the observed samples
tauMatrix = tau(1: nDiracs)';
% samples are related to locations and weights; first two terms are already
% known, solve for weights
weightsRec = locMatrix \ tauMatrix;
