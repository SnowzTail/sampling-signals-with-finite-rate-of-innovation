function [dualKernel] = dual_basis(kernel)
% Function:
%   - obtain the dual basis kernel of a biorthogonal wavelet
%
% InputArg(s):
%   - kernel: basis kernel of a biorthogonal wavelet
%
% OutputArg(s):
%   - dualKernel: the dual basis kernel of a biorthogonal wavelet
%
% Comments:
%   - calculate the inverse Gram matrix directly
%
% Author & Date: Yang (i@snowztail.com) - 15 Dec 18

% Gram matrix
gramMatrixInv = zeros(length(kernel));
for iRow = 1: length(kernel)
    for iCol = 1: length(kernel)
       gramMatrixInv(iRow, iCol) =  dot(kernel(iRow), kernel(iCol));
    end
end
% normalise the result
gramMatrixInv = gramMatrixInv / norm(gramMatrixInv);
% dual basis
dualKernel = kernel * gramMatrixInv;
end

