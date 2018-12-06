clear; close;
% max degree of polynomials
nDegMax = 3;
% kernels of finite support
lSignal = 2048;
% sampling period
tSample = 64;
% max amplitue
ampMax = 32;
% number of iterations
iter = 6;
% number of shifts
nShifts = 32;
% time of sampling points
t = 1 / tSample: 1 / tSample : nShifts;
poly = zeros(nDegMax + 1, lSignal);
coeffs = zeros(nDegMax + 1, nShifts + 1);
phi = zeros(1, lSignal);
phiSet = zeros(1 + nShifts, length(phi));
%% Daubechies
% polynomials of max degree N can be reproduced by a scaling function that
% generates wavelets with (N + 1) vanishing moments
[phiT, psiT, xval] = wavefun('dB4', iter);
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
    coeffs(iDeg + 1, :) = dot(repmat(poly(iDeg + 1, :), nShifts + 1, 1), phiSet, 2).';
end
% reproduce polynomials with coefficients
polyRep = coeffs * phiSet ./ tSample;
figure;
for iDeg = 0: nDegMax
    subplot((nDegMax + 1) / 2, 2, iDeg + 1);
    plot(t, poly(iDeg + 1, :), 'k');
    xlabel('Time');
    ylabel('Amplitude');
    hold on;
    plot(t, polyRep(iDeg + 1, :), 'r');
    xlabel('Time');
    ylabel('Amplitude');
    legend('Original', 'Reproduced');
    title(['Polynomial of degree ', num2str(iDeg)]);
end
% suptitle(['Reproduction of polynomials of maximum degree ' nDegMax ' using Daubechies wavelets with ' nDegMax + 1 ' vanishing moments']);
flag = 1;
