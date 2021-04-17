%Question 1.1
% Introduction to some simples function of image processing
clear; clc;

%Question.1.1.1
% load the picture
pic1 = imread('Bird_Q.1.1.jpg');

% RGB split - we could do with imsplit function too
redChannel = pic1(:,:,1);
greenChannel = pic1(:,:,2);
blueChannel = pic1(:,:,3);

% image and three splited channels
subplot(2,2,1);
imshow(pic1);
title('original picture','interpreter','latex');
subplot(2,2,2);
imshow(blueChannel);
title('Blue Channel','interpreter','latex');
subplot(2,2,3);
imshow(redChannel);
title('Red Channel','interpreter','latex');
subplot(2,2,4);
imshow(greenChannel);
title('green Channel','interpreter','latex');


%Question.1.1.2
figure;
sum = (blueChannel + redChannel + greenChannel);
subplot(1,2,1);
imshow(sum);
title('mean of sum of RGB channels','interpreter','latex');
subplot(1,2,2);
imshow(rgb2gray(pic1));
title('RGB to Gray using rgb2gray','interpreter','latex');

%Question.1.1.3
% make the grayscale form of the picture binary
figure;
binarypic1 = imbinarize(rgb2gray(pic1));
imshowpair(pic1,binarypic1,'montage');
title('binarized image');


