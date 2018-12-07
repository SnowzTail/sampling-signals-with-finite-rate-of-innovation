clear; close;
% max degree of polynomials
nDegMax = 3;
% kernels of finite support
lSignal = 2048;
% sampling period
tSample = 64;
% max amplitue
weightMax = 32;
% number of iterations
iter = 6;
% number of shifts
nShifts = 31;
% time of sampling points
t = 0: 1 / tSample : (lSignal - 1) / tSample;
poly = zeros(nDegMax + 1, lSignal);
coeffs = zeros(nDegMax + 1, nShifts + 1);
phi = zeros(1, lSignal);
phiSet = zeros(1 + nShifts, length(phi));
phiTilde = zeros(1, lSignal);
phiTildeSet = zeros(1 + nShifts, length(phiTilde));
%% B-splines
% polynomials of max degree N can be reproduced by a scaling function that
% generates wavelets with (N + 1) vanishing moments
[phiT, psiT, phiTildeT, psiTildeT, xval] = wavefun('bior4.4', iter);
% extend the scaling function to the signal length
phi(1: length(phiT)) = phiT;
phiTilde(1: length(phiTildeT)) = phiTildeT;
for iShift = 0: nShifts
    % obtain all shifted phi (circular shift)
    phiSet(iShift + 1, :) = circshift(phi, iShift * tSample);
    phiTildeSet(iShift + 1, :) = circshift(phiTilde, iShift * tSample);
    % obtain all shifted phi (remove tails)
%     phiSet(iShift + 1, :) = [zeros(1, iShift * tSample), phi(1: end - iShift * tSample)];
end
% dot(phiTildeSet(2,:), phiSet(1,:))
for iDeg = 0: nDegMax
    % all possible polynomials determined by time, max degree nDegMax
    poly(iDeg + 1, :) = t .^ iDeg;
    % coefficients of corresponding kernels
    coeffs(iDeg + 1, :) = dot(repmat(poly(iDeg + 1, :), nShifts + 1, 1), phiTildeSet, 2).' ./ tSample;
end
% reproduce polynomials with coefficients
polyRep = coeffs * phiSet;
figure;
for iDeg = 0: nDegMax
    subplot((nDegMax + 1) / 2, 2, iDeg + 1);
    plot(t, poly(iDeg + 1, :), 'b');
    hold on;
    plot(t, polyRep(iDeg + 1, :), 'r');
    hold on;
    plot(t, phiSet' .* coeffs(iDeg + 1, :), 'k:')
    xlabel('Time');
    ylabel('Amplitude');
    legend(sprintf('Original Polynomial t^%d', iDeg), 'Reproduced by Cubic B-spline', 'Shifted Kernels', 'location', 'northwest');
    title(['Polynomial of Degree ', num2str(iDeg)]);
end
% suptitle(['Reproduction of polynomials of maximum degree ' nDegMax ' using Cubic B-spline']);
