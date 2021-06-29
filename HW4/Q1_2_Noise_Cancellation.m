%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Q.1 - image noise cancellation
% Salt&Pepper(SP)-GaussianNoise(G)-PoissonNoise(P)-SpeckleNoise(S)
clear; clc;

% the image i`ve used is a bit large and the code takes a little time to run be patient please :) or 
% you can import a smaller photo for testing


% now each time we add these 4 noises on an image and see the result
image = imread('breakingbad.jpg');
figure;
imshow(image);
title('original image');
Noisy_Sp_image = imnoise(image,'salt & pepper');
Noisy_G_image = imnoise(image,'gaussian');
Noisy_P_image = imnoise(image,'poisson');
Noisy_S_image = imnoise(image,'speckle');
figure;
subplot(2,2,1),imshow(Noisy_Sp_image);
title('Salt&Pepper Noise added');
subplot(2,2,2),imshow(Noisy_G_image);
title('Gaussian Noise added');
subplot(2,2,3),imshow(Noisy_P_image);
title('Poisson Noise added');
subplot(2,2,4),imshow(Noisy_S_image);
title('Speckle Noise added');

% gaussian filter on 4 images
figure;
% 1 - Sp
SpGaussianfilteredImage = gaussianFilter(13,Noisy_Sp_image,1.3);
subplot(1,2,1),imshow(Noisy_Sp_image);
title('Salt&Pepper Noise added');
subplot(1,2,2),imshow(SpGaussianfilteredImage);
title('Filtered by Gaussian Filter, kernelSize = 13, sigma = 1.3');
% 2 - G
figure;
gGaussianfilteredImage = gaussianFilter(7,Noisy_G_image,0.84);
subplot(1,2,1),imshow(Noisy_G_image);
title('Gaussian Noise added');
subplot(1,2,2),imshow(gGaussianfilteredImage);
title('Filtered by Gaussian Filter, kernelSize = 7, sigma = 0.84');
% 3 - P
figure;
pGaussianfilteredImage = gaussianFilter(11,Noisy_P_image,1.1);
subplot(1,2,1),imshow(Noisy_P_image);
title('Poisson Noise added');
subplot(1,2,2),imshow(pGaussianfilteredImage);
title('Filtered by Gaussian Filter, kernelSize = 11, sigma = 1.1');
% 4 - S
figure;
sGaussianfilteredImage = gaussianFilter(15,Noisy_S_image,1.3);
subplot(1,2,1),imshow(Noisy_S_image);
title('Speckle Noise added');
subplot(1,2,2),imshow(sGaussianfilteredImage);
title('Filtered by Gaussian Filter, kernelSize = 15, sigma = 1.3');

% median filter
% 1 - Sp
SpMedianfilteredImage = medianFilter(3,Noisy_Sp_image);
figure;
subplot(1,2,1),imshow(Noisy_Sp_image);
title('Salt&Pepper Noise added');
subplot(1,2,2),imshow(SpMedianfilteredImage);
title('Filtered by Median Filter, kernelSize = 3');
% 2 - G
gMedianfilteredImage = medianFilter(5,Noisy_G_image);
figure;
subplot(1,2,1),imshow(Noisy_G_image);
title('Gaussian Noise added');
subplot(1,2,2),imshow(gMedianfilteredImage);
title('Filtered by Median Filter, kernelSize = 5');
% 3 - P
pMedianfilteredImage = medianFilter(3,Noisy_P_image);
figure;
subplot(1,2,1),imshow(Noisy_P_image);
title('Poisson Noise added');
subplot(1,2,2),imshow(pMedianfilteredImage);
title('Filtered by Median Filter, kernelSize = 3');
% 4 - S
sMedianfilteredImage = medianFilter(7,Noisy_S_image);
figure;
subplot(1,2,1),imshow(Noisy_S_image);
title('Speckle Noise added');
subplot(1,2,2),imshow(sMedianfilteredImage);
title('Filtered by Median Filter, kernelSize = 7');


% SNR of noisy images
SpNoiseSNR = snrCalculator(image,Noisy_Sp_image)
gNoiseSNR = snrCalculator(image,Noisy_G_image)
pNoiseSNR = snrCalculator(image,Noisy_P_image)
sNoiseSNR = snrCalculator(image,Noisy_S_image)

% SNR of filtered images
SpGaussianFilteredSNR = snrCalculator(image,SpGaussianfilteredImage)
gGaussianFilteredSNR = snrCalculator(image,gGaussianfilteredImage)
pGaussianFilteredSNR = snrCalculator(image,pGaussianfilteredImage)
sGaussianFilteredSNR = snrCalculator(image, sGaussianfilteredImage)

SpMedianFilteredSNR = snrCalculator(image,SpMedianfilteredImage)
gMedianFilteredSNR = snrCalculator(image,gMedianfilteredImage)
pMedianFilteredSNR = snrCalculator(image,pMedianfilteredImage)
sMedianFilteredSNR = snrCalculator(image, sMedianfilteredImage)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Q.2 - modern noise cancellation methods
clc; clear;
brainIm = imread('brain.jpg');
brainIm = im2double(rgb2gray(brainIm));
% add a gaussian noise to the image
Noisy_G_image = imnoise(brainIm,'gaussian',0,0.0025);
figure;
montage({brainIm,Noisy_G_image});
title('grayscale and noisy brain image');

% ----------------- low pass disk filters
lowBandWidthfilter = lowPassDiskFilter(62,length(Noisy_G_image));
mediumBandWidthfilter = lowPassDiskFilter(123,length(Noisy_G_image));
highBandWidthfilter = lowPassDiskFilter(246,length(Noisy_G_image));

figure;
subplot(1,3,1),imshow((lowBandWidthfilter));
title('disk lowpass filter with radius 0.2');
subplot(1,3,2),imshow((mediumBandWidthfilter));
title('disk lowpass filter with radius 0.4');
subplot(1,3,3),imshow((highBandWidthfilter));
title('disk lowpass filter with radius 0.8');

% apply these filters to the noisy image
lowBandWidthFilteredImage = abs(ifft2(ifftshift(((fftshift(fft2(Noisy_G_image)).*(lowBandWidthfilter))))));
mediumBandWidthFilteredImage = abs(ifft2(ifftshift(((fftshift(fft2(Noisy_G_image)).*(mediumBandWidthfilter))))));
highBandWidthFilteredImage = abs(ifft2(ifftshift(((fftshift(fft2(Noisy_G_image)).*(highBandWidthfilter))))));
figure;
subplot(1,3,1),imshow((lowBandWidthFilteredImage));
title('filteredImage with radius 0.2');
subplot(1,3,2),imshow((mediumBandWidthFilteredImage));
title('filteredImage with radius 0.4');
subplot(1,3,3),imshow((highBandWidthFilteredImage));
title('filteredImage with radius 0.8');

% EPI - edge preserving index
Epi_noisyIm = epiCalculator(brainIm,Noisy_G_image)
Epi_lowBWdiskFilter = epiCalculator(brainIm,lowBandWidthFilteredImage)
Epi_mediumBWdiskFilter = epiCalculator(brainIm,mediumBandWidthFilteredImage)
Epi_highBWdiskFilter = epiCalculator(brainIm,highBandWidthFilteredImage)
% SNR 
Snr_noisyIm = snrCalculator(brainIm,Noisy_G_image)
Snr_lowBWdiskFilter = snrCalculator(brainIm,lowBandWidthFilteredImage)
Snr_mediumBWdiskFilter = snrCalculator(brainIm,mediumBandWidthFilteredImage)
Snr_highBWdiskFilter = snrCalculator(brainIm,highBandWidthFilteredImage)

% -------------- gaussian filtering
gaussianFilteredImageLowVar = gaussianFilter(3,Noisy_G_image,0.4);
gaussianFilteredImageMediumVar = gaussianFilter(3,Noisy_G_image,0.8);
gaussianFilteredImageHighVar = gaussianFilter(3,Noisy_G_image,1.3);

figure;
subplot(1,3,1),imshow((gaussianFilteredImageLowVar));
title('gaussianfilteredImage with variance 0.16');
subplot(1,3,2),imshow((gaussianFilteredImageMediumVar));
title('gaussianfilteredImage with variance 0.64');
subplot(1,3,3),imshow((gaussianFilteredImageHighVar));
title('gaussianfilteredImage with variance 1.69');

% EPI - edge preserving index
Epi_gaussianFilteredLowVar = epiCalculator(brainIm,im2double(gaussianFilteredImageLowVar))
Epi_gaussianFilteredMediumVar = epiCalculator(brainIm,im2double(gaussianFilteredImageMediumVar))
Epi_gaussianFilteredHighVar = epiCalculator(brainIm,im2double(gaussianFilteredImageHighVar))
% SNR 
Snr_gaussianFilteredLowVar = snrCalculator(brainIm,im2double(gaussianFilteredImageLowVar))
Snr_gaussianFilteredMediumVar = snrCalculator(brainIm,im2double(gaussianFilteredImageMediumVar))
Snr_gaussianFilteredHighVar = snrCalculator(brainIm,im2double(gaussianFilteredImageHighVar))

% ---------------- average kernel
averageFilteredSize3 = averageFilter(3,Noisy_G_image);
averageFilteredSize5 = averageFilter(5,Noisy_G_image);
averageFilteredSize7 = averageFilter(7,Noisy_G_image);

figure;
subplot(1,3,1),imshow((averageFilteredSize3));
title('averagefilteredImage, kernel size = 3');
subplot(1,3,2),imshow((averageFilteredSize5));
title('averagefilteredImage, kernel size = 5');
subplot(1,3,3),imshow((averageFilteredSize7));
title('averagefilteredImage, kernel size = 7');

% EPI - edge preserving index
Epi_averageFilteredSize3 = epiCalculator(brainIm,(averageFilteredSize3))
Epi_averageFilteredSize5 = epiCalculator(brainIm,(averageFilteredSize5))
Epi_averageFilteredSize7 = epiCalculator(brainIm,(averageFilteredSize7))
% SNR 
Snr_averageFilteredSize3 = snrCalculator(brainIm,(averageFilteredSize3))
Snr_averageFilteredSize5 = snrCalculator(brainIm,(averageFilteredSize5))
Snr_averageFilteredSize7 = snrCalculator(brainIm,(averageFilteredSize7))

% Bilateral Filtering
NoiseSTDd = 0.05;
bilateralFilteredImage = bilateralFilter(Noisy_G_image,NoiseSTDd);
figure;
subplot(1,2,1),imshow(Noisy_G_image);
title('Noisy Image');
subplot(1,2,2),imshow(bilateralFilteredImage);
title('filtered by bilateral filter');

% EPI - edge preserving index
Epi_bilateralFiltered = epiCalculator(brainIm,(bilateralFilteredImage))
% SNR 
Snr_bilateralFiltered = snrCalculator(brainIm,(bilateralFilteredImage))

%% functions
% median filter
function outputimage = medianFilter(kernelSize,inputImage)
    % at first zeros must be padded around the image
    imageSize = size(inputImage);
    zeroPadded_Image = zeros(imageSize(1)+(kernelSize-1),imageSize(2)+(kernelSize-1),3);
    zeroPadded_Image(:,:,1) = padarray(inputImage(:,:,1),[(kernelSize-1)/2 (kernelSize-1)/2],'both');
    zeroPadded_Image(:,:,2) = padarray(inputImage(:,:,2),[(kernelSize-1)/2 (kernelSize-1)/2],'both');
    zeroPadded_Image(:,:,3) = padarray(inputImage(:,:,3),[(kernelSize-1)/2 (kernelSize-1)/2],'both');
    zeroPaded_Image = uint8(zeroPadded_Image);
    filteredImage = zeroPadded_Image;
    % moving window on 3 RGB channels
    for i=(1+(kernelSize-1)/2):(imageSize(1)+(kernelSize-1)/2)
       for j=(1+(kernelSize-1)/2):(imageSize(2)+(kernelSize-1)/2)
           % red channel
           window = zeroPadded_Image((i-(kernelSize-1)/2):(i+(kernelSize-1)/2),...
               (j-(kernelSize-1)/2):(j+(kernelSize-1)/2),1);
           vectorizedSortedWindow = sort(reshape(window,[1,kernelSize*kernelSize]));
           filteredImage(i,j,1) = vectorizedSortedWindow(1,(kernelSize*kernelSize+1)/2);
           % green channel
           window = zeroPadded_Image((i-(kernelSize-1)/2):(i+(kernelSize-1)/2),...
               (j-(kernelSize-1)/2):(j+(kernelSize-1)/2),2);
           vectorizedSortedWindow = sort(reshape(window,[1,kernelSize*kernelSize]));
           filteredImage(i,j,2) = vectorizedSortedWindow(1,(kernelSize*kernelSize+1)/2);
           % blue channel
           window = zeroPadded_Image((i-(kernelSize-1)/2):(i+(kernelSize-1)/2),...
               (j-(kernelSize-1)/2):(j+(kernelSize-1)/2),3);
           vectorizedSortedWindow = sort(reshape(window,[1,kernelSize*kernelSize]));
           filteredImage(i,j,3) = vectorizedSortedWindow(1,(kernelSize*kernelSize+1)/2);
       end
    end
    outputimage = uint8(filteredImage((1+(kernelSize-1)/2):(imageSize(1)+(kernelSize-1)/2),(1+(kernelSize-1)/2):(imageSize(2)+(kernelSize-1)/2),:));
end

% gaussian filter
function outputimage = gaussianFilter(kernelSize,inputImage,stdD)
    imageSize = size(inputImage);
    variance = stdD^2;
    % kernel
    kernel = zeros(kernelSize,kernelSize);
    for i=1:kernelSize
        for j=1:kernelSize
            kernel(i,j) = 1/(2*pi*variance)*exp(-((i-(kernelSize+1)/2)^2+(j-(kernelSize+1)/2)^2)/(2*variance));
        end
    end
    if(max(size(imageSize)) == 3) % rgb image
        filteredImage = zeros(imageSize(1)+(kernelSize-1),imageSize(2)+(kernelSize-1),3);
        filteredImage(:,:,1) = conv2(inputImage(:,:,1),kernel);
        filteredImage(:,:,2) = conv2(inputImage(:,:,2),kernel);
        filteredImage(:,:,3) = conv2(inputImage(:,:,3),kernel);
        outputimage = uint8(filteredImage((1+(kernelSize-1)/2):(imageSize(1)+(kernelSize-1)/2),(1+(kernelSize-1)/2):(imageSize(2)+(kernelSize-1)/2),:));
    end
    if(max(size(imageSize)) == 2) % grayscale image
        filteredImage = zeros(imageSize(1)+(kernelSize-1),imageSize(2)+(kernelSize-1));
        filteredImage(:,:) = conv2(inputImage,kernel);
        outputimage = filteredImage((1+(kernelSize-1)/2):(imageSize(1)+(kernelSize-1)/2),(1+(kernelSize-1)/2):(imageSize(2)+(kernelSize-1)/2));
    end
end
    

% SNR calculator
function SNR = snrCalculator(orginalImage,secondImage)
    % convert rgb to grayscale and impelement the SNR formula
    SNR = 10*log10((sum((im2gray(orginalImage)).^2,'all'))/...
        (sum((im2gray(orginalImage)-im2gray(secondImage)).^2,'all')));
end

% EPI calculator
function epi = epiCalculator(orginalImage,filteredImage)
    kernel = fspecial('laplacian',0.2);
    highP_orginalIm = conv2(orginalImage,kernel);
    highP_filteredIm = conv2(filteredImage,kernel);
    meanHighP_orginalIm = mean(highP_orginalIm(:));
    meanHighP_filteredIm = mean(highP_filteredIm(:));
    numerator = sum((highP_orginalIm-meanHighP_orginalIm).*...
        (highP_filteredIm-meanHighP_filteredIm),'all');
    denominator = sqrt((sum((highP_orginalIm-meanHighP_orginalIm).^2,'all'))*...
        (sum((highP_filteredIm-meanHighP_filteredIm).^2,'all')));
    epi = numerator/denominator;
end

% disk lowpas filter
function filter = lowPassDiskFilter(bandwidth,imageSize)
    filter = zeros(imageSize,imageSize);
    for i=1:(imageSize)
        for j=1:(imageSize)
            if ((i-(imageSize+1)/2)^2+(j-(imageSize+1)/2)^2) < (bandwidth^2)
                filter(i,j) = 1;
            end
        end
    end
end


% average filter 
function filtered = averageFilter(kernelSize,inputImage)
    imageSize = size(inputImage);
    filtered = zeros(imageSize(1)+(kernelSize-1),imageSize(2)+(kernelSize-1));
    imageSize = size(inputImage);
    kernel = ones(kernelSize,kernelSize)./kernelSize^2;
    filter = (conv2(inputImage,kernel));
    filtered = filter((1+(kernelSize-1)/2):(imageSize(1)+(kernelSize-1)/2),(1+(kernelSize-1)/2):(imageSize(2)+(kernelSize-1)/2));
end

% bilateral filter
function filtered = bilateralFilter(inputImage,noiseSTD)
    inputImage = im2double(inputImage);
    imageSize = size(inputImage);
    filtered = inputImage;
    % initialize the weight ghx
    Ghx_kernel = zeros(3,3);
    diagonalSizeIm = imageSize(1) * sqrt(2);
    hx = 0.02 * diagonalSizeIm;
    hg = 1.95 * noiseSTD;
    for i=1:3
        for j=1:3
            Ghx_kernel(i,j) = exp(-1*((i-2)^2+(j-2)^2)/(2*hx^2));
        end
    end
    % because we have medical image(MRI of brain) in our implementation, border pixels
    % are ignored because usually they don`t have any specific data
    for i=2:imageSize(1)-1
        for j=2:imageSize(1)-1
            Ghg_kernel = zeros(3,3);
            % Ghg intensity weight
            for k=-1:1
                for l=-1:1
                    Ghg_kernel(k+2,l+2) = exp(-1*((inputImage(i+k,j+l)-inputImage(i,j))^2)/(2*hg^2));
                end
            end
            filtered(i,j) = sum(inputImage((i-1):(i+1),(j-1):(j+1)).*Ghx_kernel.*Ghg_kernel,'all')/...
                sum(Ghx_kernel.*Ghg_kernel,'all');
        end
    end
end