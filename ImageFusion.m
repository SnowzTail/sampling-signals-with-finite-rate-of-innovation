function [SR_color]= ImageFusion(Tx_RGB, Ty_RGB)
% *************************************************************************
% Wavelets and Applications Course - Dr. P.L. Dragotti
% MATLAB mini-project 'Sampling Signals with Finite Rate of Innovation'
% Exercice 6
% *************************************************************************
% 
% This function fuses the set of 40 low-resolution images
% 'LR_Tiger_xx.tif' using the shifts obtained after registration and
% outputs a single super-resolved image.
%
% Input:   Tx_RGB: horizontal shifts, a 40x3 matrix
%          Ty_RGB: vertical shifts, a 40x3 matrix
%
% Tx_RGB and Ty_RGB must computed relatively to the first image in the set
% 'LR_Tiger_01.tif' and for each layer Red, Green and Blue. 
%
% note: Tx_RGB(1,:) = Ty_RGB(1,:) = (0 0 0) by definition.
%
% Author:   Loic Baboulaz
% Date:     August 2006
%
% Imperial College London
% *************************************************************************

warning off all

% Parameters
NbSensors = 40;
HR_size = 512;
LR_size = 64;
Scale = HR_size/LR_size;

% Round the input parameters if not already done
Tx_RGB = round(Tx_RGB);
Ty_RGB = round(Ty_RGB);

% load the PSF
load('2DBspline.mat','Spline2D');

HR_Samples_col = repmat(1:Scale:HR_size,LR_size,1);
HR_Samples_row = HR_Samples_col';

fprintf('\nFusing and restoring data ...')
h = waitbar(0,'Fusing and restoring data ...');
for layer = 1:3,

    waitbar(layer/3)
    SamplesOneCol = [];
    SamplesOneCol_col = [];
    SamplesOneCol_row = [];

    for k =1:NbSensors,
        % Load images
        LR_k= double(imread(sprintf('LR_Tiger_%.2d.tif',k)));
        LR_k = LR_k(:,:,layer)/255;

        SamplesOneCol = [SamplesOneCol ; reshape(LR_k ,LR_size*LR_size,1)];
        SamplesOneCol_col = [SamplesOneCol_col; reshape(HR_Samples_col - Tx_RGB(k,layer) ,LR_size*LR_size,1)];
        SamplesOneCol_row = [SamplesOneCol_row; reshape(HR_Samples_row - Ty_RGB(k,layer) ,LR_size*LR_size,1)];
    end

    % Interpolation of 2D scattered data with Matlab function griddata
    ti = 1:1:HR_size;
    [HRgrid_row HRgrid_col] = meshgrid(ti,ti);
    InterpImage(:,:,layer) = griddata(SamplesOneCol_row, SamplesOneCol_col, SamplesOneCol, HRgrid_row, HRgrid_col,'cubic').';

    % Rescale samples
    DataImage = Scale^2*InterpImage(:,:,layer);
    
    % Interpolation sometimes produces NaN values that has to be removed
    DataImage(find(isnan(DataImage)==1)) = 0;
        
    % Restore with Wiener filter
    DataImage = edgetaper(DataImage,Spline2D);
    SR = deconvwnr(DataImage,Spline2D);

    % Constraint and normalize
    SR = ((SR>=0) .*SR)/max(SR(:));

    % Store SR image
    SR_color(:,:,layer) = SR;
    
    clear SR
end
close(h);
fprintf('finished!\n');

% ----------------------- OUTPUT DATA -------------------------------------
% load ground truth image
HR_ref= double(imread('HR_Tiger_01.tif'));
HR_ref = HR_ref/255;


% load low-resolution image of reference
LR_ref= double(imread('LR_Tiger_01.tif'));
LR_ref = LR_ref/255;

% calculate the PSNR
x = SR_color(:);
y = HR_ref(:);
d = mean( mean( (x-y).^2 ) );
m = max( max(x(:)),max(y(:)) );
PSNR = 10*log10( m/d );

fprintf('\n--------------------------------')
fprintf('\n\nINFORMATION:')
fprintf('\n\nNumber of input images: %g',NbSensors)
fprintf('\n\nSize of each input image: %g x %g pixels',LR_size,LR_size)
fprintf('\n\nSize of the super-resolved image: %g x %g pixels',HR_size,HR_size)
fprintf('\n\nPSNR of the super-resolved image: %2.2f dB',PSNR)
fprintf('\n\n--------------------------------\n\n')

figure
subplot(1,2,1)
imagesc(LR_ref)
axis image
axis off
title('Input image (reference)','Fontsize',16)

subplot(1,2,2)
imagesc(SR_color)
axis image
axis off
title('Super-resolved image','Fontsize',16)
