function [phiT] = bspline(period, nOrders)
% Function:
%   - obtain the sample of bspline of a given order
%
% InputArg(s):
%   - period: sampling period
%   - nOrders: the order of bspline
%
% OutputArg(s):
%   - phiT: sample of bspline of a given order
%
% Comments:
%   - bspline of order N can reproduce polynomials of maximum degree N
%
% Author & Date: Yang (i@snowztail.com) - 15 Dec 18

% b-spline of order 0 is the box function
boxFun = ones(1, period);
prevPhiT = boxFun;
% b-spline of order N is derived by convolving the box function with itself
% for N times
if nOrders == 0
    phiT = boxFun;
else
    for iOrder = 1: nOrders
        % pad the scaling function to a power of 2 and ensure the peak in the
        % centre
        if mod(iOrder, 2)
            phiT = [0 conv(boxFun, prevPhiT)] / period;
        else
            phiT = [conv(boxFun, prevPhiT) 0] / period;
        end
        prevPhiT = phiT;
    end
end
end
