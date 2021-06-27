% Q.1 - image noise cancellation
% Salt&Pepper(SP)-GaussianNoise(G)-PoissonNoise(P)-SpeckleNoise(S)
clc; clear;
% now each time we add these 4 noises on an image and see the result
image = imread('breakingbad.jfif');

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

% median filter
filteredImage = medianFilter(3,Noisy_Sp_image);
figure;
subplot(1,2,1),imshow(Noisy_Sp_image);
title('Salt&Pepper Noise added');
subplot(1,2,2),imshow(filteredImage);
title('Filtered by Median Filter');
%% functions
% median filter
function outputimage = medianFilter(kernelSize,inputImage)
    % at first zeros must be padded around the image
    zeroPadded_Image = zeros(2400+(kernelSize-1),3600+(kernelSize-1),3);
    zeroPadded_Image(:,:,1) = padarray(inputImage(:,:,1),[(kernelSize-1)/2 (kernelSize-1)/2],'both');
    zeroPadded_Image(:,:,2) = padarray(inputImage(:,:,2),[(kernelSize-1)/2 (kernelSize-1)/2],'both');
    zeroPadded_Image(:,:,3) = padarray(inputImage(:,:,3),[(kernelSize-1)/2 (kernelSize-1)/2],'both');
    zeroPaded_Image = uint8(zeroPadded_Image);
    filteredImage = zeroPadded_Image;
    % moving window on 3 RGB channels
    for i=(1+(kernelSize-1)/2):(2400+(kernelSize-1)/2)
       for j=(1+(kernelSize-1)/2):(3600+(kernelSize-1)/2)
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
    outputimage = uint8(filteredImage((1+(kernelSize-1)/2):(2400+(kernelSize-1)/2),(1+(kernelSize-1)/2):(3600+(kernelSize-1)/2),:));
end
