clear; close;
% max degree of polynomials
degMax = 3;
% length of kernels of finite support
len = 2048;
% sampling period
period = 64;
% max amplitue
ampMax = 32;
% number of iterations
iter = log2(period);
% number of shifts
nShifts = 31;
% time of sampling points
sampPts = 0: 1 / period : (len - 1) / period;
%% B-splines
% bspline of order N can reproduce polynomials of maximum degree N
[phiT] = bspline(period, degMax);
% obtain kernel by shifting scaling function
[kernelSet] = kernel_set(len, period, nShifts, phiT);
% calculate the dual basis kernel
[dualKernel] = dual_basis(kernelSet(1, :));
% obtain dual basis kernel set by shifting
[dualKernelSet] = kernel_set(len, period, nShifts, dualKernel);
% determine polynomials and coefficients of corresponding kernels
[poly, coefs] = polynomial_coefs(len, period, nShifts, degMax, sampPts, dualKernelSet);
% reproduce polynomials with coefficients
polyRep = coefs * kernelSet;
figure;
for iDeg = 0: degMax
    subplot((degMax + 1) / 2, 2, iDeg + 1);
    plot(sampPts, poly(iDeg + 1, :), 'b');
    hold on;
    plot(sampPts, polyRep(iDeg + 1, :), 'r');
    hold on;
    plot(sampPts, kernelSet' .* coefs(iDeg + 1, :), 'k:')
    xlabel('Time');
    ylabel('Amplitude');
    legend(sprintf('Original Polynomial t^%d', iDeg), 'Reproduced by Cubic B-spline', 'Shifted Kernels', 'location', 'northwest');
    title(['Polynomial of Degree ', num2str(iDeg)]);
end
% suptitle(['Reproduction of polynomials of maximum degree ' nDegMax ' using Daubechies wavelets with ' nDegMax + 1 ' vanishing moments']);
