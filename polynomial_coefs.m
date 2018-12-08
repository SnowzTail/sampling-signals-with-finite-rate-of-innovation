function [poly, coefs] = polynomial_coefs(len, period, nShifts, degMax, sampPts, kernelSet)
% Function: 
%   - create polynomials of order up to a maximum number
%   - calculate the coefficients that the shifted kernels (i.e. components)
%   need to reconstruct the corresponding polynomials
%
% InputArg(s):
%   - len: signal length
%   - period: sampling period
%   - nShifts: number of available shifts in the signal range
%   - degMax: maximum degree of the polynomials
%   - sampPts: sampling points
%   - kernelSet: shifted kernel sets for sampling
%
% OutputArg(s):
%   - poly: value of polynomials with order 0 to degMax at the sampling
%   points
%   - coefs: the coefficients that the shifted kernels (i.e. components)
%   need to reconstruct the corresponding polynomials
%
% Comments:
%   - coefficients can be used as weights on the samples 
%
% Author & Date: Yang (i@snowztail.com) - 08 Dec 18
poly = zeros(degMax + 1, len);
coefs = zeros(degMax + 1, nShifts + 1);
for iDeg = 0: degMax
    % all possible polynomials determined by time, max degree nDegMax
    poly(iDeg + 1, :) = sampPts .^ iDeg;
    % coefficients of corresponding kernels
    coefs(iDeg + 1, :) = dot(repmat(poly(iDeg + 1, :), nShifts + 1, 1), kernelSet, 2).' ./ period;
end
end

