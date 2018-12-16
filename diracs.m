function [signal, loc, amp] = diracs(len, period, nDiracs, ampMax)
% Function: 
%   - create fixed-length sample signal with several diracs at sampling
%   points
%
% InputArg(s):
%   - len: signal length
%   - period: sampling period
%   - nDiracs: number of diracs at sampled points
%   - ampMax: maximum amplitude of samples
%
% OutputArg(s):
%   - signal: sampled signal with diracs at sampling points
%   - loc: dirac locations
%   - amp: dirac amplitudes
%
% Comments:
%   - the original signal contains only diracs
%
% Author & Date: Yang (i@snowztail.com) - 08 Dec 18
signal = zeros(1, len);
loc = sort(randperm(len, nDiracs)) / period;
amp = randperm(ampMax, nDiracs);
% loc = [8 12];
% amp = [2 3];
signal(loc * period) = amp;
end

