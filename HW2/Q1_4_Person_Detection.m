%Question.1.4 - person detection

clc;
clear;
thresh = 1; % example
peopleDetector = vision.PeopleDetector('UprightPeople_96x48','ClassificationThreshold',thresh);
I = imread('Q1_4_TAs.jpg');
[bboxes,scores] = peopleDetector(I);
I = insertObjectAnnotation(I,'rectangle',bboxes,scores);
figure;
imshow(I);
