function [Tx_RGB, Ty_RGB]= ImageRegistration
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
load('PolynomialReproduction_coef.mat', 'Coef_0_0', 'Coef_1_0', 'Coef_0_1');

% -------- include your code here -----------
% parameters
nSensors = 40;
nLayers = 3;
% noise threshold (by dichotomy trials and errors)
thrRed = 0.28;
thrGreen = 0.29;
thrBlue = 0.295;
x = zeros(nSensors, nLayers);
y = zeros(nSensors, nLayers);
Tx_RGB = zeros(nSensors, nLayers);
Ty_RGB = zeros(nSensors, nLayers);
% load low-resolution images
for iSensor = 1: nSensors
    % obtain samples
    lrSamples = double(imread(sprintf('LR_Tiger_%.2d.tif', iSensor)));
    % rescale sample values between 0 and 1
    lrSamples = lrSamples / 255;
    % colour layer: 1 -> red, 2 -> green, 3 -> blue
    for iLayer = 1: nLayers
        % samples below noise threshold are set to zero to reduce noise
        temp = lrSamples(:, :, iLayer);
        % apply red layer threshold
        if iLayer == 1
            temp(temp < thrRed) = 0;
            lrSamples(:, :, iLayer) = temp;
        % apply green layer threshold
        elseif iLayer == 2
            temp(temp < thrGreen) = 0;
            lrSamples(:, :, iLayer) = temp;
        % apply blue layer threshold
        else
            temp(temp < thrBlue) = 0;
            lrSamples(:, :, iLayer) = temp;
        end
        % compute moments
        m_00 = sum(sum(Coef_0_0 .* lrSamples(:, :, iLayer)));
        m_01 = sum(sum(Coef_0_1 .* lrSamples(:, :, iLayer)));
        m_10 = sum(sum(Coef_1_0 .* lrSamples(:, :, iLayer)));
        % barycenters
        x(iSensor, iLayer) = m_10 / m_00;
        y(iSensor, iLayer) = m_01 / m_00;
    end
    % select the first figure as reference and calculate the shift of
    % barycenters
    Tx_RGB(iSensor, :) = x(iSensor, :) - x(1, :);
    Ty_RGB(iSensor, :) = y(iSensor, :) - y(1, :);
end
end

