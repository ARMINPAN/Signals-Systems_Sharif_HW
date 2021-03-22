% all the questions are in this code and they`ve got seperated by different
% sections in the code

%%
% Ques.1.1
% Wireworld - 2d_fourColor_Cellular Automation 
% 0 - > empty , 1 - > conductor , 2 - > electron head , 3 - > electron tail
% empty : white ,conductor : green, electron head : purple , electron tail : red 
% I have defined a function which user can choose size of the table and
% and number of generations


WireWorld(20,200);

function WireWorld(size, gens)
    CurrentCells = zeros(size,size);
    nextStepCells = zeros(size,size);
    gensMemory = gens;
    drawTable(CurrentCells);
    title('Generation_0');
    
    %userInput - Initial condtion of the table
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
                    display([x,gens,y]);
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
%%
% Ques.1.2