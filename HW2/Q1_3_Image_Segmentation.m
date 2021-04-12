%Question.1.3
%%
% Question.1.3.3
% Image segmentation using K-means algorithm
% k is given by the user 
% Question 1.3.3

clc; clear;
% load the image
input_image = imread('Q1_3_flower.jpg');

% segmentation
% L = imsegkmeans(input_image,4);
% B = labeloverlay(input_image,L);
sizePic = size(input_image);
k = 4; % k is number of centroids
% generate random centroids 
x = randi(sizePic(1),1, k);
y = randi(sizePic(2),1, k);
output_image = cluster(input_image, x, y, k); % x and y are initialized cluster centers.
subplot(1,2,1);
imshow(input_image);
subplot(1,2,2);
imshow(output_image);
%% 
% Question.1.3.5
% Image segmentation using otsu algorithm


%% functions
% Question.1.3.3 - K-means algorithm
function output = cluster(input_image, x, y, k)
    % size of the photo
    picSize = size(input_image);
    % at first we have to flat the photo which means we have to change the
    % picture that is a (n, m, 3) 3D matrix to a (n*m, 3) 2D matrix
    imageData = reshape(input_image,[picSize(1)*picSize(2) , 3]);
    
    % now we get the RGBs of centroids
    centroids = zeros(k,3);
    for i = 1:k
        centroids(i,:) = input_image(x(i),y(i),:);
    end

    flag = 1;
    while flag
        % tag Label is a function which do clustering depends on the
        % RGB difference of each image`s pixel and centroids
        labels = tagLabel(imageData, centroids, k);
        
        % now we have to update the centroids
        for i = 1:k
            memory = centroids(i,:);
            sqrt(sum((memory - updateCentroids(labels, imageData, i)).^2))
            if(sqrt(sum((memory - updateCentroids(labels, imageData, i)).^2)) >= 100)
                centorids(i,:) = updateCentroids(labels, imageData, i);
                flag = 1;
            else
                flag = 0;
                break; % converged
            end
        end
        
    end
    
    imageData = updatePixels(labels, imageData, centroids);
    % updatePixels is a function which changes the RGB of the pixel
    % depends on the cluster of the pixel
    
    
    % make the image 3d matrix again from the 2d matrix
    output = reshape(imageData,[picSize(1), picSize(2) , 3]);
end

function labels = tagLabel(imageData, centroids, k)

    labels = zeros(1,length(imageData)); % a vector which saves labels
    for i = 1 : length(imageData)
        minimumDistance = 5000000; % just a variable for finding the minimum rgbDiff
        for j = 1 : k
            RGB_diff = ((cast(imageData(i,:),'double') - centroids(j,:)).^2);
            if((minimumDistance)  >= sqrt(sum(RGB_diff)))
                minimumDistance = sqrt(sum(RGB_diff));
                labels(i) = j;
            end
        end
    end
    
    
end

function updatedPic = updatePixels(labels, imageData, centroids)

    for i = 1 : length(imageData)
       imageData(i,:) = centroids(labels(i),:);
    end
    updatedPic = imageData;
end

function updated = updateCentroids(labels, input_image, j)
    label = find(labels == j);
    updated = [floor(mean(input_image(label,1))), floor(mean(input_image(label,2))), floor(mean(input_image(label,3)))];
end


% Question.1.3.5 - otsu-algorithm

