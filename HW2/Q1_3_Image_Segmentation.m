%Question.1.3
%%
% first part
% Image segmentation using K-means algorithm
% k is given by the user 
% Question 1.3.3

clc; clear;
% load the image
input_image = imread('Capture.jpg');

% segmentation
% L = imsegkmeans(input_image,4);
% B = labeloverlay(input_image,L);
sizePic = size(input_image);
k = 4;
x = randperm(sizePic(1), k);
y = randperm(sizePic(2), k);
output_image = cluster(input_image, x, y, k); % x and y are initialized cluster centers.
subplot(1,2,1);
imshow(input_image);
subplot(1,2,2);
imshow(output_image);
%% 
% second part
% Image segmentation using otsu algorithm


%% functions
% Question.1.3.3
function output = cluster(input_image, x, y, k)
    % size of the photo
    picSize = size(input_image);
    % at first we have to flat the photo which means we have to change the
    % picture that is a (n, m, 3) 3D matrix to a (n*m, 3) 2D matrix
    imageData = reshape(input_image,[picSize(1)*picSize(2) , 3]);
    
    flag = 1;
    while flag
        [x y]
        % tag Label is a function which do clustering depends on the
        % RGB difference of each image`s pixel and centroids
        labels = tagLabel(imageData, x, y, k, picSize);
        % updatePixels is a function which changes the RGB of the pixel
        % depends on the cluster of the pixel
        imageData = updatePixels(labels, imageData, input_image, x, y);

        % now have to update the centroids
        for i = 1:k
            memory = [x(i), y(i)];
            if(abs(memory - updateCentroids(labels, picSize, i)) >= 1)
                [x(i), y(i)] = updateCentroids(labels, picSize, i);
                flag = 1;
            else
                flag = 0;
            end
        end
    end
            
        % make the image 3d matrix again from the 2d matrix
        output = reshape(imageData,[picSize(1), picSize(2) , 3]);
end

function labels = tagLabel(imageData, x, y, k, imSize)

    labels = zeros(1,length(imageData)); % a vector which saves labels
    for i = 1 : length(imageData)
        minimumDistance = 50000*ones(1,3); % just a variable for finding the minimum rgbDiff
        for j = 1 : k
            RGB_diff = abs(imageData(i,:) - imageData((x(j)-1)*imSize(2)+y(j),:));
            if(sum(minimumDistance)  >= sum(RGB_diff))
                minimumDistance = RGB_diff;
                labels(i) = j;
            end
        end
    end
    
    
end

function updatedPic = updatePixels(labels, imageData, input_image, x, y)
    for i = 1 : length(imageData)
       imageData(i,:) = input_image(x(labels(i)), y(labels(i)), :);
    end
    updatedPic = imageData;
end

function [updatedx, updatedy] = updateCentroids(labels, imSize, j)
    
    label = [];
    
    label = find(labels == j);
    
    updatedx = floor(mean(mod(label,imSize(2))));
    updatedy = floor(mean(label./imSize(1)));

end
