% all the questions are in this code and they`ve got seperated by different
% sections in the code
% functions are defined at the end of the code
% it`s better if you comment questions in a section that you don`t want to
% check and just leave the one you want uncommented
%%
% Ques.1.1 and 1.2
% Wireworld - 2d_fourColor_Cellular Automaton 
% 0 - > empty , 1 - > conductor , 2 - > electron head , 3 - > electron tail
% empty : white ,conductor : green, electron head : purple , electron tail : red 
% you can choose size of the table and number of generations in the inputs
% of the function WireWorld
WireWorld(20,200);
% the parts which you may need to comment has been declared in the function
%%
% Ques.2.1 
% We want to calculate the Z transform of the x[n] given below
Z_Transform();
%% 
% Ques.2.2
%you can comment Questions you don`t need for a faster simulation
% Inverse Z transform
% Question 2.2.1 - finding zeros-poles of H1(z) and H2(z) by 
ZerosPoles();
% Question 2.2.2and and 2.2.4 - calculating partial fractions and inverse z transform
partial_fraction();
% Question 2.2.3 - calculating inverse z transform by iztrans
Inverse_Z_Transform();

%%
% Question.2.3
%you can comment Questions you don`t need for a faster simulation
% we want to find the impolse response of a LTI system which described by
% the difference equation below in different ways
% y[n] − 0.7y[n − 1] + 0.49y[n − 2] = 2x[n] − x[n − 1]
%1
ImpulseResp_ZTrans();

%2
ImpulseResp_Coefficients();

%3
ImpulseResp_Filter();

%via simulink
%%
% Question 2.4.
% making music - Karplus-Strong
% This system is described by the equation below
% y[n] = x[n] + 1/2(y[n-N] + y[n-(N+1)])
% N is the samples of x[n]
FS = 44100; % sampling frequency
% notes
notes = ['G', 'G', 'B', "D#", 'D', ...
'G', 'G', 'B', 'D', 'C', ...
'G', 'G', 'G', 'G', 'G', "G#", ...
"G#", "G#", "G#", "G#", 'G', 'G'];

% durations
note_durations = [330, 330, 490, 490, 750, ...
330, 330, 490, 490, 750, ...
330, 330, 330, 490, 490, 700, ...
330, 330, 330, 490, 490, 750];

% find notes frequencies
freq = NoteFrequency(notes);


alpha = 0.99; % a positive feedback for stability of system

% generate the notes
NoteGeneration(freq, note_durations, alpha);
%%
%you can comment Questions you don`t need for a faster simulation
% Question.3.1.1
Laplace_Transform();

% Question.3.1.2
stepResponse();

% Question.3.1.3
question3_1_3();

% Question.3.1.4
question3_1_4();
%%
% Question.3.2
% analyzing a circuit using laplace transform
Circuit();
%%
% Question 3.3
% simulink file has attached
% c1 = g/l , c2 = c/(ml) , c3 = 1/(ml^2)

%% all the functions
%Question 1
function WireWorld(size, gens)
    CurrentCells = zeros(size,size);
    nextStepCells = zeros(size,size);
    gensMemory = gens;
    drawTable(CurrentCells);
    title('Generation_0');
    
    %userInput - Initial condtion of the table
    %if you want to give initial state to the matrix from here change the 
    %CurrentCells below and comment Input function
    CurrentCells = Input(CurrentCells);
    
    
    % the periodic test case given in homework
    % if you don`t want to watch just this periodic behaviour, comment the line below and 
    % line 131 to 133 
    %, this periodic behavior is just for the one in the homework_1
    % which enter a electron head from the left after each 4 moves
    countDown = 1; % for setting the period

    
    pause(0.1)
    
    while gens > 0
        if sum(CurrentCells) == 0	% all cells are empty
			break;
        end
        

        for x=1:1:size
            for y=1:1:size
                if(CurrentCells(x, y) ~= 0)
                    %just for periodic test case in the homework - comment 2lines below for
                    %other inputs, we put a ElecHead in to the start
                    if(mod(countDown,5) == 0 && CurrentCells(1,y) == 1)
                        nextStepCells(1,y) = 2;
                    end
                    if CurrentCells(x, y) == 2 	% electron head to electron tail
                        nextStepCells(x, y) = 3;
                    end
                    if CurrentCells(x, y) == 3 	% electron tail to conductor
                        nextStepCells(x, y) = 1;
                    end
                    if CurrentCells(x, y) == 1 	% conductor to electron head or conductor
                        counter = 0;
                        for k = -1:1:1
                            for j = -1:1:1
                                % counting down number of electron head neigbors
                                if(k ~= 0 || j ~= 0)
                                    if(x == 1 && y == 1 && k > -1 && j > -1)
                                        if (CurrentCells(x+k, y+j) == 2)
                                            counter = counter + 1;
                                        end
                                    elseif(x == 1 && y == size && k > -1 && j < 1)
                                        if (CurrentCells(x+k, y+j) == 2)
                                            counter = counter + 1;
                                        end
                                    elseif(x == size && y == 1 && k < 1 &&j > -1)
                                        if (CurrentCells(x+k, y+j) == 2)
                                            counter = counter + 1;
                                        end
                                    elseif(x == size && y == size && k < 1 && j < 1)
                                        if (CurrentCells(x+k, y+j) == 2)
                                            counter = counter + 1;
                                        end  
                                    elseif(y ~= 1 && y~= size && x == 1 && k > -1)
                                        if (CurrentCells(x+k, y+j) == 2)
                                            counter = counter + 1;
                                        end  
                                    elseif(x ~= 1 && x~= size && y == 1 && j > -1)
                                        if (CurrentCells(x+k, y+j) == 2)
                                            counter = counter + 1;
                                        end  
                                     elseif(y ~= 1 && y~= size && x == size && k < 1)
                                        if (CurrentCells(x+k, y+j) == 2)
                                            counter = counter + 1;
                                        end
                                     elseif(x ~= 1 && x~= size && y == size && j < 1)
                                        if (CurrentCells(x+k, y+j) == 2)
                                            counter = counter + 1;
                                        end   
                                    elseif(x ~= size && y ~= size && x ~= 1 && y ~= 1)
                                        if (CurrentCells(x+k, y+j) == 2)
                                            counter = counter + 1;
                                        end  
                                    end
                                end
                            end
                        end
                        if(counter == 1 || counter == 2)
                            nextStepCells(x, y) = 2; % conductor to electron head
                        else
                            nextStepCells(x, y) = 1; % conductor to conductor
                        end
                    end
                end
            end
        end
        drawTable(nextStepCells);	% show next generation
		title(strcat('Generation ', num2str(gensMemory-gens+1)))
		CurrentCells =  nextStepCells;
        nextStepCells = zeros(size, size);
        % just for Periodic test case - comment it otherwise
        countDown = countDown + 1;
		pause(0.1)
		gens = gens - 1;
    end
end
function cells = Input(in)
    cells = in;
    while 1
        [x, y, button] = ginput(1);
        x = floor(x);
        y = floor(y);
        if button == 1 % left click 
            if cells(x,y) == 0  % empty to conductor
                cells(x,y) = 1; 
            elseif cells(x,y) == 1  % conductor to electron head
                cells(x,y) = 2;
            elseif cells(x,y) == 2  % electron head to electron tail
                cells(x,y) = 3;
            elseif cells(x,y) == 3  % electron tail to empty
                cells(x,y) = 0;
            end
            %now we have to update the table
            % draw table will update the table
            drawTable(cells);
            
        else % any other key to finish getting inputs
            break;
        end
    end
end
function drawTable(cells)
    % creating the table
    sizeOfAxis = length(cells);
    axis equal; grid on; cla; axis([1, sizeOfAxis+1, 1, sizeOfAxis+1]);
    set(gca, 'xtick', 1:sizeOfAxis+1, 'ytick', 1:sizeOfAxis+1);
    set(gca,'Color','#AEB1C9');
    
    for x = 1:sizeOfAxis
		for y = 1:sizeOfAxis
            if cells(x, y) == 1	% conductor
				rectangle('position', [x, y, 1, 1], 'facecolor', '#0B7942');
            elseif cells(x, y) == 2	% electron head
				rectangle('position', [x, y, 1, 1], 'facecolor', '#732661');
            elseif cells(x, y) == 3	% electron tail
				rectangle('position', [x, y, 1, 1], 'facecolor', '#A81818');
			end
		end
    end
end

%Question 2.1
function Z_Transform()
    syms z
    syms n integer
    assume(n >= 0);
    syms x(n)
    x(n) = cos(n*pi/6);

    X(z) = ztrans(x) % Z transform of the x[n] in the command window
    % we use numden function and after that solve func to find zeros and poles
    [numerator1, denominator1] = numden(X);
    zeros1 = solve(numerator1 == 0);
    poles1 = solve(denominator1 == 0);

    subplot(1,3,1);
    zplane(zeros1, poles1); % zero-pole plot for X(z)
    title('X(z)','interpreter','latex');


    % part 2 - Question 2.1
    X2(z) = X(2*z) % X(2z)
    % we use numden function and after that solve func to find zeros and poles
    [numerator2, denominator2] = numden(X2);
    zeros2 = solve(numerator2 == 0);
    poles2 = solve(denominator2 == 0);

    subplot(1,3,2);
    zplane(zeros2, poles2); % zero-pole plot for X(2z)
    title('X(2z)','interpreter','latex');

    % inverse z transform for X2
    xinverse1(n) = iztrans(X2);
    % change syms function to a numerical function
    xinv1_func = double(xinverse1(0:20));

    % part 3 - Question 2.1
    X3 = X(z^3)% X(z^3)
    % we use numden function and after that solve func to find zeros and poles
    [numerator3, denominator3] = numden(X3);
    zeros3 = solve(numerator3 == 0);
    poles3 = solve(denominator3 == 0);


    subplot(1,3,3);
    zplane(double(zeros3), poles3); % zero-pole plot for X(z^3)
    title('X(z power 3)','interpreter','latex');


    % inverse z transforms for X3
    xinverse2(n) = iztrans(X3);
    % change syms function to a numerical function
    xinv2_func = double(xinverse2(0:40));

    % inverse z transforms` plot
    % Time Domain functions for X(2z) and X(z^3)
    figure;
    subplot(2,1,1);
    stem([0:20],xinv1_func,'LineWidth',2,'color','b');
    title('Time domain for X(2z)','interpreter','latex');
    grid on;
    subplot(2,1,2);
    stem([0:40],xinv2_func,'LineWidth',2,'color','b');
    title('Time domain for X(z power 3)','interpreter','latex');
    grid on;
end

%Question 2.2
%Question 2.2.1 - finding zeros-poles of H1(z) and H2(z)
function ZerosPoles()
    figure;
    syms z
    H1 = (1-z^-1)/(1-z^(-1)+0.5*z^(-2));
    H2 = z^(-1)/(2-3^(1/2)*z^(-1)+0.5*z^(-2));
    [numerator1, denominator1] = numden(H1);
    zeros1 = solve(numerator1 == 0);
    poles1 = solve(denominator1 == 0);

    [numerator2, denominator2] = numden(H2);
    zeros2 = solve(numerator2 == 0);
    poles2 = solve(denominator2 == 0);
    subplot(2,1,1);
    zplane(double(zeros1),double(poles1)); % zero-pole plot for H1(z)
    title('H1(z)','interpreter','latex');
    subplot(2,1,2);
    zplane(double(zeros2),double(poles2)); % zero-pole plot for H2(z)
    title('H2(z)','interpreter','latex');
end

%Question 2.2.2 and 2.2.4 - calculating partial fractions and z inverse 2.2.2-> assume our
%system is right sided (casual) and  2.2.4-> our system is left sided and anti-casual
function partial_fraction()
    % by residuez function
    figure;
    z = tf('z');
    H1 = (1-z^-1)/(1-z^(-1)+0.5*z^(-2));
    H2 = z^(-1)/(2-3^(1/2)*z^(-1)+0.5*z^(-2));

    % at first we have to find coefficients of our transfer function using tfdata
    [numerator1, denominator1] = tfdata(H1);
    numerator1 = cell2mat(numerator1);
    denominator1 = cell2mat(denominator1);
    [numerator2, denominator2] = tfdata(H2);
    numerator2 = cell2mat(numerator2);
    denominator2 = cell2mat(denominator2);

    % remove zeros in numerator1 and denominator1 except one zero that is b0
    % in H2
    %H1 coefficients
    numerator1 = numerator1(numerator1~=0);
    denominator1 = denominator1(denominator1~=0);
    %H2 coefficients
    numerator2 = numerator2(numerator2~=0);
    denominator2 = denominator2(denominator2~=0);

    %partial fractions
    [ro1,po1,ko1] = residuez(numerator1, denominator1)
    [ro2,po2,ko2] = residuez([0,numerator2], denominator2)

    % now we have the valuse of poles and residuez so the time domain function
    % are just I have written in stem
    % now z inverse calculating - assume the inverses are right sided
    n = 0:20;
    subplot(1,2,1);
    stem([0:20],((0.5000 + 0.5000i)*(0.5000 + 0.5000i).^n + (0.5000 - 0.5000i)*(0.5000 - 0.5000i).^n),'LineWidth',2,'color','b');
    title('Casual Time Domain of H1(z) by partial fraction','interpreter','latex');
    grid on;
    subplot(1,2,2);
    stem([0:20],((-i)*(0.4330 + 0.2500i).^n + (i)*(0.4330 - 0.2500i).^n),'LineWidth',2,'color','b');
    title('Casual Time Domain of H2(z) by partial fraction','interpreter','latex');
    grid on;

    % now assume our system is anti-casual
    figure;
    n = -30:-1;
    subplot(1,2,1);
    stem([-30:-1],(-(0.5000 + 0.5000i)*(0.5000 + 0.5000i).^(n) - (0.5000 - 0.5000i)*(0.5000 - 0.5000i).^(n)),'LineWidth',2,'color','b');
    title('Anti-Casual Time Domain of H1(z) by partial fraction','interpreter','latex');
    grid on;
    subplot(1,2,2);
    stem([-30:-1],-((1i)*(0.4330 + 0.2500i).^(n) - (1i)*(0.4330 - 0.2500i).^(n)),'LineWidth',2,'color','b');
    title('Anti-Casual Time Domain of H2(z) by partial fraction','interpreter','latex');
    grid on;
end

%Question 2.2.3 - calculating inverse z transform by iztrans
function Inverse_Z_Transform()
    figure;
    syms n integer
    syms z
    H1 = (1-z^-1)/(1-z^(-1)+0.5*z^(-2));
    H2 = z^(-1)/(2-3^(1/2)*z^(-1)+0.5*z^(-2));
    %inverses
    xinverse1(n) = iztrans(H1);
    xinverse2(n) = iztrans(H2);

    %plotting the inverses
    subplot(1,2,1);
    stem([0:20],double(xinverse1(0:20)),'LineWidth',2,'color','b')
    title('Time Domain of H1(z) by iztrans','interpreter','latex');
    grid on;
    subplot(1,2,2);
    stem([0:20],double(xinverse2(0:20)),'LineWidth',2,'color','b')
    title('Time Domain H2(z) by iztrans','interpreter','latex');
    grid on;
end

% Question 2.3
% we want to find the impolse response of a LTI system which described by
% the difference equation below
% y[n] − 0.7y[n − 1] + 0.49y[n − 2] = 2x[n] − x[n − 1]
% 2.3.1 -> Method 1 - by Z Transform
function ImpulseResp_ZTrans()
    % if we get z transform of both sides, Transfer function H(z)=Y(z)/X(z) is:
    figure;

    z = tf('z');
    H = (2-z^(-1))/(1-0.7*z^(-1)+0.49*z^(-2));

    % at first we have to find coefficients of our transfer function using tfdata
    [numerator, denominator] = tfdata(H);
    numerator = cell2mat(numerator);
    denominator = cell2mat(denominator);
    % remove zeros in numerator and denominator
    %H coefficients
    numerator = numerator(numerator~=0);
    denominator = denominator(denominator~=0);

    %partial fractions
    [ro,po,ko] = residuez(numerator, denominator)

    %now we have residuez and the poles and so we have partial fraction of H(z)
    % and so we have the impulse response of the system
    n = 0:20;
    stem(n,(1.0000 + 0.2474i)*((0.3500 + 0.6062i).^n) + (1.0000 - 0.2474i)*((0.3500 - 0.6062i).^n),'LineWidth',2,'color','b');
    title('Impulse response by partial fraction','interpreter','latex');
    grid on;
end

% 2.3.2 -> Method 2 - h[n] = sigma (alpha_i * p_i^n)u[n]-finding alpha_i s with
% help of the intial conditions 
% if we calculate on paper h[0] = 2 and h[1] = 0.4
% we have the poles from part 2.3.1 which are 0.3500 +- 0.6062i
function ImpulseResp_Coefficients()
    figure;
    syms n integer;
    assume(n >= 0);
    syms alpha_1 alpha_2
    h(n) = alpha_1*(0.3500 + 0.6062i).^n + alpha_2*(0.3500 - 0.6062i).^n;
    coefs = solve(h(0) == 2, h(1) == 0.4);

    alpha1 = coefs.alpha_1
    alpha2 = coefs.alpha_2

    %impulse response
    m = 0:20;
    stem(m,alpha1*(0.3500 + 0.6062i).^m + alpha2*(0.3500 - 0.6062i).^m,'LineWidth',2,'color','b');
    title('Impulse response Question2.3.2','interpreter','latex');
    grid on;
end

% 2.3.3 -> Method 3 - using filter function to find impulse response
function ImpulseResp_Filter()
    figure;

    z = tf('z');
    H = (2-z^(-1))/(1-0.7*z^(-1)+0.49*z^(-2));

    % at first we have to find coefficients of our transfer function using tfdata
    [numerator, denominator] = tfdata(H);
    numerator = cell2mat(numerator);
    denominator = cell2mat(denominator);
    % remove zeros in numerator and denominator
    %H coefficients
    numerator = numerator(numerator~=0);
    denominator = denominator(denominator~=0);

    x = zeros(1,20);
    x(1) = 1;

    h = filter(numerator,denominator,x);
    n = 0:19;
    stem(n,h,'LineWidth',2,'color','b');
    title('Impulse response using filter function','interpreter','latex');
    grid on;
end

% 2.3.4 -> Method 4 - using simulink and designing the system
% file has been attached in the same directory


% Question 2.4
% a function which founds the frequency of notes
function freqss = NoteFrequency(note)
    len = length(note);
    freqs = zeros(1,len);
    for i=1:len
        if(note(1,i) == 'G')
            freqs(i) = 783.99;
        end
        if(note(1,i) == "G#")
            freqs(i) = 830.61;
        end
        if(note(1,i) == 'B')
            freqs(i) = 493.88;
        end
        if(note(1,i) == 'E')
            freqs(i) = 659.25;
        end
        if(note(1,i) == 'D')
            freqs(i) = 587.33;
        end
        if(note(1,i) == "D#")
            freqs(i) = 622.25;
        end
        if(note(1,i) == 'A')
            freqs(i) = 440;
        end
        if(note(1,i) == "A#")
            freqs(i) = 466.16;
        end
        if(note(1,i) == 'F')
            freqs(i) = 698.46;
        end
        if(note(i) == "F#")
            freqs(i) = 739.99;
        end
        if(note(1,i) == 'C')
            freqs(i) = 523.25;
        end
        if(note(i) == "C#")
            freqs(1,i) = 554.37;
        end
    end
    freqss = freqs;
end

% generate the notes
function NoteGeneration(freq, duration, alpha)
    noteNum = length(freq);
    FS = 44100;
    %generating notes
    for i=1:noteNum
        N = floor(FS/freq(i)); % number of samples for each note
        % *80 is just for making the music slower and more peaceful :)
        outputLength = (N*duration(i)*80);
        x = zeros(1,N);
        x(1,:) = randi([-1,1],1,N);
        x=[x zeros(1,outputLength-N)];

        y = zeros(1,outputLength); %output note
        y(1:N) = x(1,1:N);
        % transfer function
        a=[1 zeros(1,N-1) -alpha/2 -alpha/2];  %denominator 
        b=1;  %numerator 
        y = filter(b,a,x); %output
        sound(y,FS) % use to play the note
    end
end

% Question 3.1.1 - Laplace Transform
function Laplace_Transform()
    % plot poles of the systems below
    figure;

    syms s
    H1 = 1/(s^3 + 40*s^2 + 10*s + 500);
    H2 = 1/(s^4 + 12.5*s^3 + 10*s^2 + 1);

    % now we have to find the poles which are the roots of denominators
    [numerator1, denominator1] = numden(H1);
    zeros1 = [];
    poles1 = solve(denominator1 == 0);

    [numerator2, denominator2] = numden(H2);
    zeros2 = [];
    poles2 = solve(denominator2 == 0);

    % plot the poles
    subplot(1,2,1);
    zplane(zeros1,double(poles1));
    title('H1(s)','interpreter','latex');
    subplot(1,2,2);
    zplane(zeros2,double(poles2));
    title('H2(s)','interpreter','latex');
end

% Question 3.1.2 - find the response of the systems when input is unit step
function stepResponse()
    figure;

    syms s;
    syms n integer;
    assume(n >= 0);

    H1 = 1/(s^3 + 40*s^2 + 10*s + 500);
    H2 = 1/(s^4 + 12.5*s^3 + 10*s^2 + 1);

    % use ilaplace to find unit step response by inverse laplace transform
    h1_inv(n) = ilaplace(H1/s,s,n);
    h2_inv(n) = ilaplace(H2/s,s,n);

    h1 = double(h1_inv(1:0.1:100));
    h2 = double(h2_inv(1:0.1:100));

    subplot(2,1,1);
    plot((1:0.1:100),h1,'LineWidth',2,'color','b');
    title('Unit Step Response of y1(t)','interpreter','latex');
    grid on;
    subplot(2,1,2);
    plot((1:0.1:100),h2,'LineWidth',2,'color','b');
    title('Unit Step Response of y2(t)','interpreter','latex');
    grid on;
end

% Question 3.1.3 
function question3_1_3()
    figure;

    syms s;
    syms n integer;
    assume(n >= 0);

    a = [4, 5, 6];
    H1 = (2*s+1)/(s^2 + a(1)*s + 7);
    H2 = (2*s+1)/(s^2 + a(2)*s + 7);
    H3 = (2*s+1)/(s^2 + a(3)*s + 7);

    %impulse responses of the 3 systems
    h1_inv(n) = ilaplace(H1,s,n);
    h2_inv(n) = ilaplace(H2,s,n);
    h3_inv(n) = ilaplace(H3,s,n);

    h1 = double(h1_inv(0:0.01:10));
    h2 = double(h2_inv(0:0.01:10));
    h3 = double(h3_inv(0:0.01:10));

    %plot impulse responses

    plot((0:0.01:10),h1);
    hold on;
    plot((0:0.01:10),h2);
    hold on;
    plot((0:0.01:10),h3);
    title('Impulse Response of y(t) for different valuse of a','interpreter','latex');
    legend('a = 4','a = 5','a = 6');
    grid on;

    %find unit step responses
    figure;
    %unit step response responses of the 3 systems
    h1_inv(n) = ilaplace(H1/s,s,n);
    h2_inv(n) = ilaplace(H2/s,s,n);
    h3_inv(n) = ilaplace(H3/s,s,n);

    h1 = double(h1_inv(0:0.01:10));
    h2 = double(h2_inv(0:0.01:10));
    h3 = double(h3_inv(0:0.01:10));

    %plot unit step responses

    plot((0:0.01:10),h1);
    hold on;
    plot((0:0.01:10),h2);
    hold on;
    plot((0:0.01:10),h3);
    title('Unit Step Response of y(t) for different valuse of a','interpreter','latex');
    legend('a = 4','a = 5','a = 6');
    grid on;
end


% Question 3.1.4
function question3_1_4()
    figure;
    
    syms s;
    syms n integer;
    assume(n >= 0);
    H1 = tf([1 1],[1 3 4]);
    H2 = tf([1 1],[1 3 4 0]);
    H3 = tf([1 1],[1 3 4 0 0]);
   
    % use feedback function to create a feedback 
    H1_new = feedback(H1,1);
    H2_new = feedback(H2,1);
    H3_new = feedback(H3,1);
    
    % change this transfer func`s to syms functions to use in ilaplace
    [Num,Den] = tfdata(H1_new);
    H1_neww = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s);
    [Num,Den] = tfdata(H2_new);
    H2_neww = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s);
    [Num,Den] = tfdata(H3_new);
    H3_neww = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s);
    
    % unit step responses
    h1_inv(n) = ilaplace(H1_neww/s,s,n);
    h2_inv(n) = ilaplace(H2_neww/s,s,n);
    h3_inv(n) = ilaplace(H3_neww/s,s,n);
    
    h1 = double(h1_inv(0:0.1:50));
    h2 = double(h2_inv(0:0.1:50));
    h3 = double(h3_inv(0:0.1:50));
    
    %plot unit step responses
    
    subplot(1,2,1);
    plot((0:0.1:50),h1);
    hold on;
    plot((0:0.1:50),h2);
    hold on;
    plot((0:0.1:50),h3);
    title('Unit Step Response of y1(t)/y2(t)/y3(t)','interpreter','latex');
    legend('y1(t)','y2(t)','y3(t)');
    grid on;
    
    %plot ramp resposes
    
    subplot(1,2,2);
    n = 0:0.1:50;  % the time vector
    input = n;   % assuming ramp with slope=1
    y = lsim(H1_new, input, n);
    plot(n,y);
    hold on;
    y = lsim(H2_new, input, n);
    plot(n,y);
    hold on;
    y = lsim(H3_new, input, n);
    plot(n,y);
    title('Ramp Response of y1(t)/y2(t)/y3(t)','interpreter','latex');
    legend('y1(t)','y2(t)','y3(t)');
    grid on;
end

% Question 3.2
% the differential equation which describes the Q of capacitor:
% 300u(t) = 0.06 d^2(Vc)/d^2t + 0.32 dVc/dt + Vc
function Circuit()
    figure;
    
    syms s;
    syms n integer;
    assume (n>=0);
    V = 300/(s*(0.06*s^2 + 0.32*s + 1));
    V_inverse(n) = ilaplace(V,s,n);
    Q = 0.02 * V_inverse;
    Ic = 0.02 * diff(V_inverse);
    Q = double(Q(0:0.01:5));
    Ic = double(Ic(0:0.01:5));
    
    subplot(1,2,1); 
    plot((0:0.01:5),Q,'LineWidth',2,'color','b');
    title('Charge of the Capacitor','interpreter','latex');
    grid on;
    subplot(1,2,2);
    plot((0:0.01:5),Ic,'LineWidth',2,'color','b');
    title('Current of the Circuit','interpreter','latex');
    grid on;
    
end

% Question 3.3
% simulink file has attached