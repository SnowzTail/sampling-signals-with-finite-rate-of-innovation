clear; close;
%% Annihilating filter
load('tau.mat');
% max degree of polynomials
nDegMax = length(tau) - 1;
% number of pulses should be known
nDiracs = 2;
tauMatrixLeft = zeros(nDegMax - nDiracs + 1, nDiracs);
tauMatrixRight = zeros(nDegMax - nDiracs + 1, 1);
% Yule-Walker system
for iDirac = 1: nDiracs
   tauMatrixLeft(iDirac, :) = flip(tau(iDirac: iDirac + nDiracs - 1));
   tauMatrixRight(iDirac) = -tau(iDirac + nDiracs);
end
% solve equations for filter coefficients
filterCoeffs = [1; tauMatrixLeft \ tauMatrixRight];
% roots of z-transform of the filter corresponds to the pulse locations:
% H(z) = (1-t0z^(-1))(1-t1z^(-1)) = 1-(t0+t1)z^(-1)+t0t1z(-2)
% h(0)=1; h(1)=t0+t1; h(2)=t0t1
func = @(tDelta) ([tDelta(1) * tDelta(2) - filterCoeffs(3); -(tDelta(1) + tDelta(2)) - filterCoeffs(2)]);
% solve to get pulse locations
tDelta = sort(fsolve(func, [1 1]));
% Vandermonde system
timeMatrix = fliplr(vander(tDelta))';
% weighted sum of the observed samples
tauMatrix = tau(1: nDiracs)';
% samples are related to locations and weights; first two terms are already
% known, solve for weights
ampMatrix = timeMatrix \ tauMatrix;
