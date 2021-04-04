%Question.1.2
% introduction to an edge detection algorithm for digital images called
% Soble - Feldman which is a first order method
 
clear; clc;

%Question.1.1.2
% load the image
pic1 = imread('NY_Q.1.2.jpg');
%imshow(pic1);
% define the kernels
kernel_1 = [-1 0 1; -2 0 2; -1 0 1];
kernel_2 = [-1 -2 -1; 0 0 0; 1 2 1];

% use rgb2gray two find intensity of each pixel 
% comment rgb2gray if you`re photo is already grayscale and just leave pic1
Gx = conv2(kernel_1, rgb2gray(pic1));
Gy = conv2(kernel_2, rgb2gray(pic1));

G = sqrt(Gx.^2 + Gy.^2);

% normalize and print the image
imshow(G/max(G(:)));
