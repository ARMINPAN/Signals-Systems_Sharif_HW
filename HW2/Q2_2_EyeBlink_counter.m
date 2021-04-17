% Question.2.2
% count number of eye blinks in the given data

clc;
clear;
load('eye.mat');
eyeData = E; % change the variable name to a better one :)

% we want to analyze the data with a moving window with the length of
% 500 samples
% total length of our data is 300000 samples

windowLength = 500;
blinkNum = 0; % a variable to save number of blinks
blinkTime = []; % a vector to save times of blinkings
for i=1:(length(eyeData)/windowLength)
    [maximum, index] = max(eyeData((i-1)*windowLength+1:i*windowLength));
    if(maximum > 2.5) % condition of blinking is if the maximum data in a window is more than 2.5
        blinkNum = blinkNum + 1;
        blinkTime = [blinkTime, index+(i-1)*500];
    end
end
size(blinkTime)
plot(eyeData);
hold on;
line([blinkTime;blinkTime],[0;5].*ones(size(blinkTime)),'MarkerFaceColor','#A2142F','Color','#A2142F','Marker','.');
