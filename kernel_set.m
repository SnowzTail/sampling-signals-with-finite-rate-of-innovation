function [kernelSet] = kernel_set(len, period, nShifts, phiT)
% Function: 
%   - obtain the shifted kernel sets for sampling at corresponding points
%
% InputArg(s):
%   - len: signal length
%   - period: sampling period
%   - nShifts: number of available shifts in the signal range
%   - phiT: scaling function with finite support (shorter than signals)
%
% OutputArg(s):
%   - kernelSet: shifted kernel sets for sampling
%
% Comments:
%   - the scaling function should be determined in advance

% Author & Date: Yang (i@snowztail.com) - 08 Dec 18
kernel = zeros(1, len);
kernelSet = zeros(1 + nShifts, len);
% extend the scaling function to the signal length
kernel(1: length(phiT)) = phiT;
for iShift = 0: nShifts
    % obtain all shifted phi (circular shift)
%     kernelSet(iShift + 1, :) = circshift(kernel, iShift * period);
    % obtain all shifted phi (remove tails)
    kernelSet(iShift + 1, :) = [zeros(1, iShift * period), kernel(1: end - iShift * period)];
end
end

