clear; close;
% number of diracs
nDiracs = 2;
% max degree of polynomials
degMax = ceil(2 * nDiracs - 1);
% kernels of finite support
len = 2048;
% sampling period
period = 64;
% max amplitue
ampMax = 32;
% number of shifts
nShifts = 31;
% number of iterations
iter = log2(period);
% time of sampling points
% t = 1 / tSample: 1 / tSample : lSignal / tSample;
sampPts = 0: 1 / period : (len - 1) / period;
%% Daubechies
% polynomials of max degree N can be reproduced by a scaling function that
% generates wavelets with (N + 1) vanishing moments
[phiT, ~, ~] = wavefun('dB4', iter);
% obtain kernel by shifting scaling function
[kernelSet] = kernel_set(len, period, nShifts, phiT);
% determine polynomials and coefficients of corresponding kernels
[poly, coefs] = polynomial_coefs(len, period, nShifts, degMax, sampPts, kernelSet);
%% Recover signal from samples
% samples are known
load ('samples.mat');
samples = y_sampled;
tau = zeros(1, degMax + 1);
for iDeg = 0: degMax
    tau(1, iDeg + 1) = dot(coefs(iDeg + 1, :), samples);
end
tauMatrixLeft = zeros(degMax - nDiracs + 1, nDiracs);
tauMatrixRight = zeros(degMax - nDiracs + 1, 1);
% Yule-Walker system
for iDirac = 1: degMax - nDiracs + 1
   tauMatrixLeft(iDirac, :) = flip(tau(iDirac: iDirac + nDiracs - 1));
   tauMatrixRight(iDirac) = -tau(iDirac + nDiracs);
end
% solve equations for filter coefficients
filterCoefs = [1; tauMatrixLeft \ tauMatrixRight];
% roots of the filter corresponds to the pulse locations
locEst = sort(zero(tf(filterCoefs',1)))';
% Vandermonde system
locMatrix = fliplr(vander(locEst))';
% weighted sum of the observed samples
tauMatrix = tau(1: nDiracs)';
% samples are related to locations and weights; first two terms are already
% known, solve for weights
ampEst = (locMatrix \ tauMatrix)';
%% Signal plot
stem(locEst, ampEst, 'r--x');
xlabel('Time');
ylabel('Amplitude');
legend('Estimated Diracs');
title('Reconstruction of Dirac Signal with Kernels Reproducing Polynomials');
% disp(['Exact location:      ' num2str(loc)]);
disp(['Estimated location:  ' num2str(locEst)]);
% disp(['Exact amplitude:     ' num2str(amp)]);
disp(['Estimated amplitude: ' num2str(ampEst)]);
