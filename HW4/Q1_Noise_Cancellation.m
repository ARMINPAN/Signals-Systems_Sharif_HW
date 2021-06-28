% Q.1 - image noise cancellation
% Salt&Pepper(SP)-GaussianNoise(G)-PoissonNoise(P)-SpeckleNoise(S)

% the image i`ve used is a bit large and the code takes a little time to run be patient please :) or 
% you can import a smaller photo for testing


% now each time we add these 4 noises on an image and see the result
image = imread('breakingbad.jpg');

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
SpGaussianfilteredImage = gaussianFilter(7,Noisy_Sp_image,0.84);
subplot(1,2,1),imshow(Noisy_Sp_image);
title('Salt&Pepper Noise added');
subplot(1,2,2),imshow(SpGaussianfilteredImage);
title('Filtered by Gaussian Filter');
% 2 - G
figure;
gGaussianfilteredImage = gaussianFilter(7,Noisy_G_image,0.84);
subplot(1,2,1),imshow(Noisy_G_image);
title('Gaussian Noise added');
subplot(1,2,2),imshow(gGaussianfilteredImage);
title('Filtered by Gaussian Filter');
% 3 - P
figure;
pGaussianfilteredImage = gaussianFilter(7,Noisy_P_image,0.84);
subplot(1,2,1),imshow(Noisy_P_image);
title('Poisson Noise added');
subplot(1,2,2),imshow(pGaussianfilteredImage);
title('Filtered by Gaussian Filter');
% 4 - S
figure;
sGaussianfilteredImage = gaussianFilter(7,Noisy_S_image,0.84);
subplot(1,2,1),imshow(Noisy_S_image);
title('Speckle Noise added');
subplot(1,2,2),imshow(sGaussianfilteredImage);
title('Filtered by Gaussian Filter');

% median filter
% 1 - Sp
SpMedianfilteredImage = medianFilter(3,Noisy_Sp_image);
figure;
subplot(1,2,1),imshow(Noisy_Sp_image);
title('Salt&Pepper Noise added');
subplot(1,2,2),imshow(SpMedianfilteredImage);
title('Filtered by Median Filter');
% 2 - G
gMedianfilteredImage = medianFilter(3,Noisy_G_image);
figure;
subplot(1,2,1),imshow(Noisy_G_image);
title('Gaussian Noise added');
subplot(1,2,2),imshow(gMedianfilteredImage);
title('Filtered by Median Filter');
% 3 - P
pMedianfilteredImage = medianFilter(3,Noisy_P_image);
figure;
subplot(1,2,1),imshow(Noisy_P_image);
title('Poisson Noise added');
subplot(1,2,2),imshow(pMedianfilteredImage);
title('Filtered by Median Filter');
% 4 - S
sMedianfilteredImage = medianFilter(3,Noisy_S_image);
figure;
subplot(1,2,1),imshow(Noisy_S_image);
title('Speckle Noise added');
subplot(1,2,2),imshow(sMedianfilteredImage);
title('Filtered by Median Filter');


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
% Q.2 - modern noise cancellation methods



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
    % zero padding
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
          filteredImage(i,j,1) = sum(kernel.*zeroPadded_Image((i-(kernelSize-1)/2):(i+(kernelSize-1)/2),...
               (j-(kernelSize-1)/2):(j+(kernelSize-1)/2),1),'all');
          % green channel
          filteredImage(i,j,2) = sum(kernel.*zeroPadded_Image((i-(kernelSize-1)/2):(i+(kernelSize-1)/2),...
               (j-(kernelSize-1)/2):(j+(kernelSize-1)/2),2),'all');
           % blue channel
          filteredImage(i,j,3) = sum(kernel.*zeroPadded_Image((i-(kernelSize-1)/2):(i+(kernelSize-1)/2),...
               (j-(kernelSize-1)/2):(j+(kernelSize-1)/2),3),'all');
       end
    end
    outputimage = uint8(filteredImage((1+(kernelSize-1)/2):(imageSize(1)+(kernelSize-1)/2),(1+(kernelSize-1)/2):(imageSize(2)+(kernelSize-1)/2),:));
end
    

% SNR calculator
function SNR = snrCalculator(orginalImage,secondImage)
    % convert rgb to grayscale and impelement the SNR formula
    SNR = 10*log10((sum((im2gray(orginalImage)).^2,'all'))/...
        (sum((im2gray(orginalImage)-im2gray(secondImage)).^2,'all')));
end


