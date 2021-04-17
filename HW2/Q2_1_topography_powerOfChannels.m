% Topography and power of channels
% Question.2.1

clc; 
clear;

load('data.mat');
time = [5; 40; 20; 40; 20; 40; 20; 40; 20; 40; 20; 40; 5];
% sampling frequency is 250Hz
% the first and the last 5 seconds is resting state and we remove them
data = data(1251:length(data)-1250,:);
time = time(2:12);

% now we normalize the data 
data = zscore(data);

% we have 5 rests with the length of 20secs
rests = zeros(100*250,19);
summ = zeros(20*250,19);

for i=1:5
    rests(((i-1)*20*250)+1:((i*20*250)),:) = data((sum(time(1:(2*i-1)))*250+1):(sum(time(1:(2*i-1)))*250+20*250),:);
    summ = summ + rests(((i-1)*20*250)+1:((i*20*250)),:);
end

% mean of rest times for all channels  with length of 20secs\5000
rests = summ/5; % rest signal


% we have 6 task with the length of 40secs
tasks = zeros(240*250,19);
summ = zeros(40*250,19);

for i=1:6
    tasks(((i-1)*40*250)+1:((i*40*250)),:) = data(((sum(time(1:(2*i-1)))*250+1)-40*250):(sum(time(1:((2*i-1))))*250),:);
    summ = summ + tasks(((i-1)*40*250)+1:((i*40*250)),:);
end

% mean of rest times for all channels  with length of 40secs\10000
tasks = summ/5; % task signal


% power of each electrod for rest and task signal
Prest = 1/length(rests)*sum(rests(:,:).^2,1);
Ptask = 1/length(rests)*sum(tasks(:,:).^2,1);

% now we plot the tophography map
channel_title = {'FP1','FP2','F7','F3','FZ','F4','F8','T7','C3','CZ','C4','T8','P7','P3','PZ','P4','P8','O1','O2'};

figure;
% use min for a better plot
plot_topography(channel_title,min(Prest,0.25),1,'10-20',1,1,1000); % rest state
title('topography of resting state');
figure;
plot_topography(channel_title,min(Ptask,0.48),1,'10-20',1,1,1000); % task state
title('topography of task state');