function [Tx_RGB Ty_RGB]= ImageRegistration
% *************************************************************************
% Wavelets and Applications Course - Dr. P.L. Dragotti
% MATLAB mini-project 'Sampling Signals with Finite Rate of Innovation'
% Exercice 6
% *************************************************************************
% 
% FOR STUDENTS
%
% This function registers the set of 40 low-resolution images
% 'LR_Tiger_xx.tif' and returns the shifts for each image and each layer
% Red, Green and Blue. The shifts are calculated relatively to the first
% image 'LR_Tiger_01.tif'. Each low-resolution image is 64 x64 pixels.
%
%
% OUTPUT:   Tx_RGB: horizontal shifts, a 40x3 matrix
%           Ty_RGB: vertical shifts, a 40x3 matrix
%
% NOTE: _Tx_RGB(1,:) = Ty_RGB(1,:) = (0 0 0) by definition.
%       _Tx_RGB(20,2) is the horizontal shift of the Green layer of the
%       20th image relatively to the Green layer of the firs image.
%
%
% OUTLINE OF THE ALGORITHM:
%
% 1.The first step is to compute the continuous moments m_00, m_01 and m_10
% of each low-resolution image using the .mat file called:
% PolynomialReproduction_coef.mat. This file contains three matrices
% 'Coef_0_0', 'Coef_1_0' and 'Coef_0_1' used to calculate the continuous
% moments.
%
% 2.The second step consists in calculating the barycenters of the Red,
% Green and Blue layers of the low-resolution images.
%
% 3.By computing the difference between the barycenters of corresponding 
% layers between two images, the horizontal and vertical shifts can be 
% retrieved for each layer.
%
%
% Author:   Loic Baboulaz
% Date:     August 2006
%
% Imperial College London
% *************************************************************************


% Load the coefficients for polynomial reproduction
load('PolynomialReproduction_coef.mat','Coef_0_0','Coef_1_0','Coef_0_1');

% -------- include your code here -----------
NbSensors = 40;
% for layer = 1:3
for k=1:NbSensors
    % Load images
    LR_k= double(imread(sprintf('LR_Tiger_%.2d.tif',k)));
    LR_k = LR_k(:,:,:)/255;
    for layer = 1:3
    m00 = sum(sum(Coef_0_0 .* LR_k(:, :, layer)));
    m01 = sum(sum(Coef_0_1 .* LR_k(:, :, layer)));
    m10 = sum(sum(Coef_1_0 .* LR_k(:, :, layer)));
    x(layer, k) = m10/m00;
    y(layer, k) = m01/m00;
    end
end
for k=1:NbSensors
    LR_k= double(imread(sprintf('LR_Tiger_%.2d.tif',k)));
    LR_k = LR_k(:,:,:)/255;
    Tx_RGB = x(:, k) - x(:, 1);
    Ty_RGB = y(:, k) - y(:, 1);
end

end
% end

