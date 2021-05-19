% audio signal processing
% question.1
%% part 1 - DTMF
%  import the filters in workspace at first to run code
clc; clear;
load('BP_697.mat');
load('BP_770.mat');
load('BP_852.mat');
load('BP_941.mat');
load('BP_1209.mat'); 
load('BP_1336.mat'); 
load('BP_1477.mat'); 
load('BP_1633.mat'); 


clc;
% for other recorded signals change the input of audioread:
% DialedSequence_NoNoise.wav
% DialedSequence_SNR00dB.wav
% DialedSequence_SNR10dB.wav
% DialedSequence_SNR20dB.wav
% DialedSequence_SNR30dB.wav

[audio, Fs1] = audioread('DialedSequence_NoNoise.wav');

% apply the filters
bp697 = conv(BP_697, audio);
bp770 = conv(BP_770, audio);
bp852 = conv(BP_852, audio);
bp941 = conv(BP_941, audio);
bp1209 = conv(BP_1209, audio);
bp1336 = conv(BP_1336, audio);
bp1477 = conv(BP_1477, audio);
bp1633 = conv(BP_1633, audio);

Filtered_signal = bp697+bp770+bp852+bp941+bp1209+bp1336+bp1477+bp1633;

subplot(2,1,1);
plot(audio);
title('Original Signal','interpreter','latex');
xlabel('Sample','interpreter','latex');
ylabel('Amplitude','interpreter','latex');
subplot(2,1,2);
plot(Filtered_signal);
title('Filtered Signal','interpreter','latex');
xlabel('Sample','interpreter','latex');
ylabel('Amplitude','interpreter','latex');

% detect whether a DTMF is presented at a particullar time or not by
% setting a threshold C 
% if the signal is more than the threshold in any time, we relize that a
% key is pressed
% by looking at the smoothed plots we realize that C = 0.15 is a good
% threshold 

bp697_pressed = find(bp697 > 0.15);
bp770_pressed = find(bp770 > 0.15);
bp852_pressed = find(bp852 > 0.15);
bp941_pressed = find(bp941 > 0.15);
bp1209_pressed = find(bp1209 > 0.15);
bp1336_pressed = find(bp1336 > 0.15);
bp1477_pressed = find(bp1477 > 0.15);
bp1663_pressed = find(bp1633 > 0.15);

% find intersection of the sets
one_pressed = min(intersect(bp697_pressed, bp1209_pressed));
two_pressed = min(intersect(bp697_pressed, bp1336_pressed));
three_pressed = min(intersect(bp697_pressed, bp1477_pressed));
A_pressed = min(intersect(bp697_pressed, bp1663_pressed));

four_pressed = min(intersect(bp770_pressed, bp1209_pressed));
five_pressed = min(intersect(bp770_pressed, bp1336_pressed));
six_pressed = min(intersect(bp770_pressed, bp1477_pressed));
B_pressed = min(intersect(bp770_pressed, bp1663_pressed));

seven_pressed = min(intersect(bp852_pressed, bp1209_pressed));
eight_pressed = min(intersect(bp852_pressed, bp1336_pressed));
nine_pressed = min(intersect(bp852_pressed, bp1477_pressed));
C_pressed = min(intersect(bp852_pressed, bp1663_pressed));

star_pressed = min(intersect(bp941_pressed, bp1209_pressed));
zero_pressed = min(intersect(bp941_pressed, bp1336_pressed));
hashtag_pressed = min(intersect(bp941_pressed, bp1477_pressed));
D_pressed = min(intersect(bp941_pressed, bp1663_pressed));

% order of pressed keys
Order = sort([one_pressed, two_pressed, three_pressed, A_pressed, ...
    four_pressed, five_pressed, six_pressed, B_pressed, seven_pressed, ...
    eight_pressed, nine_pressed, C_pressed, star_pressed, zero_pressed, ...
    hashtag_pressed, D_pressed]);

% print the decoded data
decoded = [];
for i=1:length(Order)
    if Order(i) == one_pressed 
        fprintf("1");
        decoded = [decoded, '1'];
    elseif Order(i) == two_pressed 
        fprintf("2");
        decoded = [decoded, '2'];
    elseif Order(i) == three_pressed 
        fprintf("3");
        decoded = [decoded, '3'];
    elseif Order(i) == A_pressed 
        fprintf("A");
        decoded = [decoded, 'A'];
    elseif Order(i) == four_pressed 
        fprintf("4");
        decoded = [decoded, '4'];
    elseif Order(i) == five_pressed 
        fprintf("5");
        decoded = [decoded, '5'];
    elseif Order(i) == six_pressed 
        fprintf("6");
        decoded = [decoded, '6'];
    elseif Order(i) == B_pressed 
        fprintf("B");
        decoded = [decoded, 'B'];
    elseif Order(i) == seven_pressed 
        fprintf("7");
        decoded = [decoded, '7'];
    elseif Order(i) == eight_pressed 
        fprintf("8");
        decoded = [decoded, '8'];
    elseif Order(i) == nine_pressed 
        fprintf("9");
        decoded = [decoded, '9'];
    elseif Order(i) == C_pressed 
        fprintf("C");
        decoded = [decoded, 'C'];
    elseif Order(i) == star_pressed 
        fprintf("*");
        decoded = [decoded, '*'];
    elseif Order(i) == zero_pressed 
        fprintf("0");
        decoded = [decoded, '0'];
    elseif Order(i) == hashtag_pressed 
        fprintf("#");
        decoded = [decoded, '#'];
    elseif Order(i) == D_pressed 
        fprintf("D");
        decoded = [decoded, 'D'];
    end
    
end

% part 2 - a GUI for decoded keys using attached - visualized dial

ShowDecodedKeys(decoded);

function ShowDecodedKeys(input)
keys = {}; % output image

for i=1:length(input)
    if(input(i) == '1')
        keys{end+1} =  imread('1.bmp');
    elseif(input(i) == '2')
        keys{end+1} = imread('2.bmp');
    elseif(input(i) == '3')
        keys{end+1} = imread('3.bmp');
    elseif(input(i) == 'A')
        keys{end+1} = imread('A.bmp');    
    elseif(input(i) == '4')
        keys{end+1} =  imread('4.bmp');   
    elseif(input(i) == '5')
        keys{end+1} = imread('5.bmp');   
    elseif(input(i) == '6')
       keys{end+1} = imread('6.bmp');   
    elseif(input(i) == 'B')
        keys{end+1} =  imread('B.bmp');   
    elseif(input(i) == '7')
        keys{end+1} = imread('7.bmp');       
    elseif(input(i) == '8')
        keys{end+1} = imread('8.bmp');  
    elseif(input(i) == '9')
        keys{end+1} = imread('9.bmp');   
    elseif(input(i) == 'C')
        keys{end+1} = imread('C.bmp');  
    elseif(input(i) == '*')
        keys{end+1} = imread('asterisco.bmp');  
    elseif(input(i) == '0')
        keys{end+1} = imread('0.bmp'); 
    elseif(input(i) == '#')
        keys{end+1} = imread('cancelletto.bmp');
    elseif(input(i) == 'D')
        keys{end+1} = imread('D.bmp');
    end
end

figure;
montage(keys,'Size',[1 length(input)]);
title('Visualized Dials','interpreter','latex');
end




