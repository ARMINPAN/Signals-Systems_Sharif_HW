% all the questions are in this code and they`ve got seperated by different
% sections in the code
% functions are defined at the end of the code
%%
% Ques.1.1 and 1.2
% 1.1
% Wireworld - 2d_fourColor_Cellular Automation 
% 0 - > empty , 1 - > conductor , 2 - > electron head , 3 - > electron tail
% empty : white ,conductor : green, electron head : purple , electron tail : red 
% I have defined a function which user can choose size of the table and
% and number of generations
WireWorld(20,200);

%%
% Ques.2.1
% We want to calculate the Z transform of the x[n] given below
Z_Transform();
%%
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
    
    pause(0.01)
    
    while gens > 0
        if sum(CurrentCells) == 0	% all cells are empty
			break;
        end
        for x=1:1:size
            for y=1:1:size
                if(CurrentCells(x, y) ~= 0)
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
                                    elseif(x == 1 && k > -1)
                                        if (CurrentCells(x+k, y+j) == 2)
                                            counter = counter + 1;
                                        end  
                                    elseif(y == 1 && j > -1)
                                        if (CurrentCells(x+k, y+j) == 2)
                                            counter = counter + 1;
                                        end  
                                     elseif(x == size && k < 1)
                                        if (CurrentCells(x+k, y+j) == 2)
                                            counter = counter + 1;
                                        end
                                     elseif(y == size && j < 1)
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
		pause(0.01)
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

X(z) = ztrans(x); % Z transform of the x[n] in the command window
% we use numden function and after that solve func to find zeros and poles
[numerator1, denominator1] = numden(X);
zeros1 = solve(numerator1 == 0);
poles1 = solve(denominator1 == 0);

subplot(1,3,1);
zplane(zeros1, poles1); % zero-pole plot for X(z)
title('X(z)');


% part 2 - Question 2.1
X2(z) = X(2*z); % X(2z)
% we use numden function and after that solve func to find zeros and poles
[numerator2, denominator2] = numden(X2);
zeros2 = solve(numerator2 == 0);
poles2 = solve(denominator2 == 0);

subplot(1,3,2);
zplane(zeros2, poles2); % zero-pole plot for X(2z)
title('X(2z)');

% inverse z transform for X2
xinverse1(n) = iztrans(X2);


% part 3 - Question 2.1
X3 = X(z^3);% X(z^3)
% we use numden function and after that solve func to find zeros and poles
[numerator3, denominator3] = numden(X3);
zeros3 = solve(numerator3 == 0);
poles3 = solve(denominator3 == 0);


subplot(1,3,3);
zplane(zeros3(3:6), poles3); % zero-pole plot for X(z^3)
title('X(z^3)');
%I have just used zeros(3:6) because ine the numerator of X3 there is a z^3
%and z^3 = 0 would put three zeros in out vector and this cause errors if
%we use zeros3(1:6) while plotting zeros

% inverse z transforms for X3
xinverse2(n) = iztrans(X3);

% inverse z transforms` plot
figure;
subplot(2,1,1);
ezplot(n,xinverse1);
subplot(2,1,2);
ezplot(n,xinverse2);

 
end

%Qu
