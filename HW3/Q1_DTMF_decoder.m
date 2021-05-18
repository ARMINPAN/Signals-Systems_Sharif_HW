% audio signal processing
% question.1

%% part 1
clc;
%1. no noise audio
[audioNoNoise, Fs1] = audioread('DialedSequence_NoNoise.wav');

% apply the filters
bp697 = conv(BP_697.Numerator, audioNoNoise);
bp770 = conv(BP_770.Numerator, audioNoNoise);
bp852 = conv(BP_852.Numerator, audioNoNoise);
bp941 = conv(BP_941.Numerator, audioNoNoise);
bp1209 = conv(BP_1209.Numerator, audioNoNoise);
bp1336 = conv(BP_1336.Numerator, audioNoNoise);
bp1477 = conv(BP_1477.Numerator, audioNoNoise);
bp1633 = conv(BP_1633.Numerator, audioNoNoise);

Filtered_signal = bp697+bp770+bp852+bp941+bp1209+bp1336+bp1477+bp1633;

plot(Filtered_signal);

% detect whether a DTMF is presented at a particullar time or not by
% setting a threshold C 
% if the signal is more than the threshold in any time, we relize that a
% key is pressed
% by looking at the smoothed plots we realize that C = 0.25 is a good
% threshold for no noise audio


bp697_pressed = find(bp697 > 0.3);
bp770_pressed = find(bp770 > 0.3);
bp852_pressed = find(bp852 > 0.3);
bp941_pressed = find(bp941 > 0.3);
bp1209_pressed = find(bp1209 > 0.3);
bp1336_pressed = find(bp1336 > 0.3);
bp1477_pressed = find(bp1477 > 0.3);
bp1663_pressed = find(bp1633 > 0.3);

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

Order = sort([one_pressed, two_pressed, three_pressed, A_pressed, ...
    four_pressed, five_pressed, six_pressed, B_pressed, seven_pressed, ...
    eight_pressed, nine_pressed, C_pressed, star_pressed, zero_pressed, ...
    hashtag_pressed, D_pressed]);

% print the decoded data
for i=1:length(Order)
    if Order(i) == one_pressed 
        fprintf("1");
    elseif Order(i) == two_pressed 
        fprintf("2");
    elseif Order(i) == three_pressed 
        fprintf("3");
    elseif Order(i) == A_pressed 
        fprintf("A");
    elseif Order(i) == four_pressed 
        fprintf("4");
    elseif Order(i) == five_pressed 
        fprintf("5");
    elseif Order(i) == six_pressed 
        fprintf("6");
    elseif Order(i) == B_pressed 
        fprintf("B");
    elseif Order(i) == seven_pressed 
        fprintf("7");
    elseif Order(i) == eight_pressed 
        fprintf("8");
    elseif Order(i) == nine_pressed 
        fprintf("9");
    elseif Order(i) == C_pressed 
        fprintf("C");
    elseif Order(i) == star_pressed 
        fprintf("*");
    elseif Order(i) == zero_pressed 
        fprintf("0");
    elseif Order(i) == hashtag_pressed 
        fprintf("#");
    elseif Order(i) == D_pressed 
        fprintf("D");
    end
    
end