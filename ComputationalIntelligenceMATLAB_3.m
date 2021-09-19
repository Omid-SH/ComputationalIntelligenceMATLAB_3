%% Data Visualization
clear; clc
x = load('DataNew.mat');
x = x.DataNew;
scatter(x(1,:), x(2,:), 'filled');

%% Genetic clustering

% first method (not good result)
init_vector = [ones(1,200), 2.*ones(1,200), 3.*ones(1,200), 4.*ones(1,200), 5.*ones(1,200)];
init_vector = init_vector(randperm(1000));
options = optimoptions(options,'InitialPopulationMatrix', init_vector);
options = optimoptions(options,'FitnessScalingFcn', @fitscalingprop);
options = optimoptions(options,'SelectionFcn', @selectionroulette);
options = optimoptions(options,'CrossoverFcn', @crossovertwopoint);
options = optimoptions(options,'Display', 'off');
options = optimoptions(options,'PlotFcn', {  @gaplotbestf @gaplotscorediversity @gaplotscores });

[result, fval] = ga(@genetic_maximum_five_fitness,1000,[],[],[],[],ones(1,1000),5.*ones(1,1000),[],1:1000,options);
scatter(x(1,:),x(2,:),15,result, 'filled');
title(['fval maximum 5 cluster GA = ' num2str(fval)]);


%% maximum 5 clusters

options = optimoptions('ga');
options = optimoptions(options,'FitnessScalingFcn', @fitscalingprop);
options = optimoptions(options,'SelectionFcn', @selectionroulette);
options = optimoptions(options,'CrossoverFcn', @crossoverarithmetic);
options = optimoptions(options,'MutationFcn', {  @mutationgaussian 1  0.7 });
options = optimoptions(options,'Display', 'off');
options = optimoptions(options,'PlotFcn', {  @gaplotbestf @gaplotscorediversity @gaplotscores });

[result, fval] = ga(@genetic_maximum_five_fitness_1,10,[],[],[],[],-40.*ones(1,10),40.*ones(1,10),[],[],options);

% Show result
c = ones(1,1000);
for i = 1:1000 
    minn = norm([result(1), result(2)] - [x(1,i), x(2,i)]);
    for j = 2:5
        if norm([result(2*j-1), result(2*j)] - [x(1,i), x(2,i)]) < minn
            minn = norm([result(2*j-1), result(2*j)] - [x(1,i), x(2,i)]);
            c(i) = j;
        end
    end
end
figure()
scatter(x(1,:),x(2,:),15,c, 'filled');
title(['fval maximum 5 cluster GA = ' num2str(fval)]);

    
%% exactly 5 clusters

options = optimoptions('ga');
options = optimoptions(options,'FitnessScalingFcn', @fitscalingprop);
options = optimoptions(options,'SelectionFcn', @selectionroulette);
options = optimoptions(options,'CrossoverFcn', @crossoverarithmetic);
options = optimoptions(options,'MutationFcn', {  @mutationgaussian 1  0.4 });
options = optimoptions(options,'Display', 'off');
options = optimoptions(options,'PlotFcn', {  @gaplotbestf @gaplotscorediversity @gaplotscores });

[result, fval] = ga(@genetic_exactly_five_fitness_1,10,[],[],[],[],-40.*ones(1,10),40.*ones(1,10),[],options);

% Show result
c = ones(1,1000);
for i = 1:1000 
    minn = norm([result(1), result(2)] - [x(1,i), x(2,i)]);
    for j = 2:5
        if norm([result(2*j-1), result(2*j)] - [x(1,i), x(2,i)]) < minn
            minn = norm([result(2*j-1), result(2*j)] - [x(1,i), x(2,i)]);
            c(i) = j;
        end
    end
end
figure()
scatter(x(1,:),x(2,:),15,c, 'filled');
title(['fval exactly 5 cluster GA = ' num2str(fval)]);

    
%% PSO clustering

%% maximum 5 clusters
options = optimoptions('particleswarm');
options = optimoptions(options,'PlotFcn', {  @pswplotbestf });
[result, fval] = particleswarm(@genetic_maximum_five_fitness_1, 10, -40.*ones(1,10), 40.*ones(1,10), options);

% Show result
c = ones(1,1000);
for i = 1:1000 
    minn = norm([result(1), result(2)] - [x(1,i), x(2,i)]);
    for j = 2:5
        if norm([result(2*j-1), result(2*j)] - [x(1,i), x(2,i)]) < minn
            minn = norm([result(2*j-1), result(2*j)] - [x(1,i), x(2,i)]);
            c(i) = j;
        end
    end
end
figure()
scatter(x(1,:),x(2,:),15,c, 'filled');
title(['fval maximum 5 cluster PSO = ' num2str(fval)]);

%% exactly 5 clusters

options = optimoptions('particleswarm');
options = optimoptions(options,'PlotFcn', {  @pswplotbestf });

[result, fval] = particleswarm(@genetic_exactly_five_fitness_1, 10, -40.*ones(1,10), 40.*ones(1,10), options);

% Show result
c = ones(1,1000);
for i = 1:1000 
    minn = norm([result(1), result(2)] - [x(1,i), x(2,i)]);
    for j = 2:5
        if norm([result(2*j-1), result(2*j)] - [x(1,i), x(2,i)]) < minn
            minn = norm([result(2*j-1), result(2*j)] - [x(1,i), x(2,i)]);
            c(i) = j;
        end
    end
end
figure()
scatter(x(1,:),x(2,:),15,c, 'filled');
title(['fval exactly 5 cluster PSO = ' num2str(fval)]);

%% ACO clustering

%% maximum 5 clusters
nIteration = 1000;
nAnt = 100;
fEvapapiration = 0.9;
fG = 1.*ones(5, 1000);

c_best = zeros(1, 1000);
cost_best = inf;
for i = 1:nIteration
    
    costs = zeros(1, nAnt); 
    c = zeros(nAnt, 1000); % cluster number each ant select for each node

    for ant=1:nAnt
        s = 0; % sum of
        for j = 1:1000
            s = sum(fG(:,j));
            r = s.*rand();
            p = 0;
            for k = 1:5
                p = p + fG(k,j);
                if( r <= p )
                    c(ant, j) = k;
                    break;
                end
            end
        end

        % Prevent pheromone distribution by sorting and renaming labels
        c(ant, :) = reorder(c(ant, :));
        
        % Calculate ant's answer cost and update fG  
        cost = genetic_maximum_five_fitness(c(ant, :));
        costs(ant) = cost;
        
        % Update best answer
        if(cost < cost_best)
            c_best = c(ant, :);
            cost_best = cost;
        end
        
    end

    for ant=1:nAnt
        fG = fG .* fEvapapiration; % Evapapiration
        f = 10^13 ./ (costs(ant)^3);

        % Update pheromone
        for j = 1:1000
            fG(c(ant, j), j) = fG(c(ant, j), j) + f;
        end
    end
end

figure()
scatter(x(1,:),x(2,:),15,c_best(1,:), 'filled');
title(['fval maximum 5 cluster ACO = ' num2str(cost_best)]);

%% exactly 5 clusters

nIteration = 1000;
nAnt = 100;
fEvapapiration = 0.9;
fG = 1.*ones(5, 1000);

c_best = zeros(1, 1000);
cost_best = inf;

for i = 1:nIteration
    
    costs = zeros(1, nAnt); 
    c = zeros(nAnt, 1000); % cluster number each ant select for each node

    for ant=1:nAnt
        s = 0; % sum of
        for j = 1:1000
            s = sum(fG(:,j));
            r = s.*rand();
            p = 0;
            for k = 1:5
                p = p + fG(k,j);
                if( r <= p )
                    c(ant, j) = k;
                    break;
                end
            end
        end

        % Prevent pheromone distribution by sorting and renaming labels
        c(ant, :) = reorder(c(ant, :));
        
        % Calculate ant's answer cost and update fG    
        cost = genetic_exactly_five_fitness(c(ant, :));
        costs(ant) = cost;
        
        % Update best answer
        if(cost < cost_best)
            c_best = c(ant, :);
            cost_best = cost;
        end
    end

    for ant=1:nAnt
        fG = fG .* fEvapapiration; % Evapapiration
        f = 10^13 ./ (costs(ant)^3);

        % Update pheromone
        for j = 1:1000
            fG(c(ant, j), j) = fG(c(ant, j), j) + f;
        end
    end
end

figure()
scatter(x(1,:),x(2,:),15,c(1,:), 'filled');
title(['fval exactly 5 cluster ACO = ' num2str(genetic_maximum_five_fitness(c(1,:)))]);

%% ACO Extra (Code from Github repo: https://github.com/madvn/ACO)

clc;
clear all;
close all;
 
%% load data 
 
data = load('DataNew.mat');
data = data.DataNew;
data = data';

dl = length(data);
%plot(data(:,1),data(:,2),'*');
 
%% init
colonySize = 50;
xAxis = 1:(colonySize);
yAxis = 1:(colonySize);
 
% parameters from data
mu = mean(data);
covariance = cov(data);
 
%visibility
s = 8;
%iterations
N = 2500;
alpha = .35;
kp = 0.2;
kd = 0.05;
phimax = 5;
phimin = 0.1;
 
food = 1:dl;
antCount = 30;
foodMat = zeros(colonySize);
antMat = zeros(colonySize);
foodLoc = zeros(dl,2);
antLoc = zeros(antCount,2);
pheromoneMat = phimin*ones(colonySize);
pheromoneIncerement = .5;
pheromoneDecrement = .99;
antMouths = zeros(antCount,1);
 
figure;
hold on;
grid on;
 
%% create ants and drop data points(food) randomly in grid
%  %  disp('Let there be an ants');
 
for i=1:antCount
    exit = 0;
    while(exit == 0)
        xPos = ceil(rand*colonySize);
        yPos = ceil(rand*colonySize);
        flag = 0;
        if(antMat((xPos),(yPos))== 0 && xPos ~=0 && yPos ~= 0)
            antMat((xPos),(yPos)) = i;
            antLoc(i,1) = xPos;
            antLoc(i,2) = yPos;
            exit = 1;
        end
    end
end
 
% %  %  disp('done dropping ants');
plot(antLoc(:,1),antLoc(:,2),'r*');
for i=1:dl
    exit = 0;
    while(exit == 0)
        xPos = ceil(rand*colonySize);
        yPos = ceil(rand*colonySize);
        if(foodMat(ceil(xPos),ceil(yPos))== 0 && xPos ~=0 && yPos ~= 0)
            foodMat(ceil(xPos),ceil(yPos)) = i;
            foodLoc(i,1) = xPos;
            foodLoc(i,2) = yPos;
            exit = 1;
        end
    end
end
plot(foodLoc(:,1),foodLoc(:,2),'go');
pickMat = zeros(N,antCount);
 
%% choose direction to move
 
for iter = 1:N
    for ant = 1:antCount
        % move
        xPos = antLoc(ant,1);
        yPos = antLoc(ant,2);
 
        % compute possible paths
        dirOptions = [xPos yPos+1; xPos+1 yPos+1; xPos+1 yPos; xPos+1 yPos-1; xPos yPos-1; xPos-1 yPos-1; xPos-1 yPos; xPos-1 yPos+1];
        %  %  disp('getting possible directions')
        dirLim = ones(1,8);
        % find available directions with limiting conditions
        if (xPos+1 > colonySize)
            dirLim([2 3 4]) = 0;
        end
        if (xPos-1 <= 0)
            dirLim([6 7 8]) = 0;
        end
        if (yPos+1 > colonySize)
            dirLim([1 2 8]) = 0;
        end
        if(yPos-1 <= 0)
            dirLim([4 5 6]) = 0;
        end
        %  %  disp('level1')
        % shd check antLoc for positions - because if there is another ant
        % in a prospective postion this ant cannot move there
        for t = 1:8
            if(dirLim(t) == 1)
                if(find(antLoc(:,1) == dirOptions(t,1)))
                    if(find(antLoc(:,2) == dirOptions(t,2)))
                        dirLim(t) = 0;
                    end
                end
            end
        end
        %  %  disp('level2');
        %if ant has food in mouth it cant move to a position that has food
        if(antMouths(ant,1) ~= 0)
            for t = 1:8
                if(dirLim(t) == 1)
                    if(foodMat(dirOptions(t,1),dirOptions(t,2)) ~= 0)
                        dirLim(t) = 0;
                    end
                end
            end
        end
        % creating a matrix that holds the indeces of the possible
        % directions of movement
        %  %  disp('levlel3');
        validDirIndeces = zeros(length(find(dirLim)),1);
        len = 0;
        for t = 1:8
            if(dirLim(t) == 1)
                validDirIndeces(len+1,1) = t;
                len = len+1;
            end
        end
        %  %  disp('final level');
        % select a direction
        if(find(validDirIndeces))
 
            % compute prob along each possible direction
            dirProbs = zeros(length(validDirIndeces),1);
            %  %  disp('dirProbs');
            for dir=1:length(validDirIndeces)
%                 dirProbs(dir,1) = pheromoneMat(dirOptions(validDirIndeces(dir,1),1),dirOptions(validDirIndeces(dir,1),2))*rand*rand;
                  dirProbs(dir,1) = rand;
            end
 
            maxProbIndex = find(dirProbs == max(dirProbs));
            %  %  disp('maxProbIndex');
            % set co-ordinates to move
            move2 = zeros(1,2);
            move2 = dirOptions(validDirIndeces(maxProbIndex(1),1),:);
            %  %  disp(move2);
% move
            antMat(move2(1,:)) = ant;
            antLoc(ant,:)=move2(1,:);
            if(antMouths(ant,1) ~= 0)
                foodMat(xPos,yPos) = 0;
                foodMat(move2(1,1),move2(1,2))=antMouths(ant,1);
                foodLoc(antMouths(ant,1),:) = move2(1,:);
            end
%             pheromoneMat(move2(1,1),move2(1,2)) = pheromoneMat(move2(1,1),move2(1,2)) + pheromoneIncerement;
%             if(pheromoneMat(move2(1,1),move2(1,2)) > phimax)
%                 pheromoneMat(move2(1,1),move2(1,2)) = phimax;
%             end
            %  %  disp('incremented pheromone');
            nxPos = move2(1,1);
            nyPos = move2(1,2);
            % pick up/drop
            foodOptions = [nxPos nyPos+1; nxPos+1 nyPos+1; nxPos+1 nyPos; nxPos+1 nyPos-1; nxPos nyPos-1; nxPos-1 nyPos-1; nxPos-1 nyPos; nxPos-1 nyPos+1];
            % find valid food directions out the above options
            foodDirLim = ones(1,8);
            % find available directions with limiting conditions
            if (nxPos+1 > colonySize)
                foodDirLim([2 3 4]) = 0;
            end
            if (nxPos-1 <= 0)
                foodDirLim([6 7 8]) = 0;
            end
            if (nyPos+1 > colonySize)
                foodDirLim([1 2 8]) = 0;
            end
            if(nyPos-1 <= 0)
                foodDirLim([4 5 6]) = 0;
            end
% pick up/drop food
            % check if mouth is empty
            %  %  disp('created valid foodDirs');
            if(antMouths(ant,1) ~= 0)
                % no?
                %  %  disp('mouth not empty');
                % see if it can be dropped
                % calculate distance between food in mouth and surrounding
                % food if any
                len = 0;
                loc = zeros(8,1);
                distance(1)=0;
                f = 0;
                for t=1:8
                    if(foodDirLim(1,t) ~=0)
                        if(foodMat(foodOptions(t,1),foodOptions(t,2)) ~= 0)
                            loc(t,1) = 1;
                            distance(len+1) = euclidianDistance(data(foodMat(foodOptions(t,1),foodOptions(t,2)),:),data(antMouths(ant,1),:));
                            f = f + (1-(distance(len+1)/alpha));
                            len = len+1;
                            %  %  disp('distCal');
                        end
                    end
                end
                f = f/(s*s);
                ppick = (1/(pheromoneMat(nxPos,nyPos)*f))*(kp/(kp+f))^2;
                pdrop = (pheromoneMat(nxPos,nyPos)*f)*(f/(kd+f))^2;
%                 if(mean(distance) < distThresh|| sum(distance)> distMax)
                if(pdrop > ppick)
                    % drop it
                    foodMat(nxPos,nyPos) = antMouths(ant,1);
                    foodLoc(antMouths(ant,1),:) = [nxPos nyPos];
                    antMouths(ant,1) = 0;
                    pheromoneMat(nxPos,nxPos) = pheromoneMat(nxPos,nxPos) + pheromoneIncerement*pheromoneMat(nxPos,nxPos);
                    if(pheromoneMat(nxPos,nxPos) > phimax)
                        pheromoneMat(nxPos,nxPos) = phimax;
                    end
                     %  disp('dropped');
                end
            else
                % yes?
                % if food is present see if food can be picked up'
                %                   %  disp('mouth empty');
                if(foodMat(nxPos,nyPos) ~= 0)
                    len = 0;
                    loc = zeros(8,1);
                    distance(1)=0;
%                      %  disp('food located');
                    f = 0;
                    for t=1:8
                        if(foodDirLim(1,t) ~=0)
                            if(foodMat(foodOptions(t,1),foodOptions(t,2)) ~= 0)
                                loc(t,1) = 1;
                                distance(len+1) = euclidianDistance(data(foodMat(foodOptions(t,1),foodOptions(t,2)),:),data(foodMat(nxPos,nyPos),:));
                                f = f + (1-(distance(len+1)/alpha));
                                len = len+1;
                            end
                        end
                    end
                    if(length(distance) >= 1)
%                          %  disp(mean(distance));
                        f = f/(s*s);
                        ppick = (1/(pheromoneMat(nxPos,nyPos)*f))*(kp/(kp+f))^2;
                        pdrop = (pheromoneMat(nxPos,nyPos)*f)*(f/(kd+f))^2;
                        if(ppick > pdrop)
                            % pick it up
                            antMouths(ant,1) = foodMat(nxPos,nyPos);
                            foodMat(nxPos,nyPos) = 0;
%                              %  disp('pickedUp');
                            pheromoneMat(nxPos,nxPos) = pheromoneDecrement*pheromoneMat(nxPos,nxPos);
                            if(pheromoneMat(nxPos,nxPos) < phimin)
                                pheromoneMat(nxPos,nxPos) = phimin;
                            end
                        end
                    end
                end
            end
        end
% update map for location of food and ants
        pause(.01);
        cla;
        plot(foodLoc(:,1),foodLoc(:,2),'go');
        plot(antLoc(:,1),antLoc(:,2),'r*');
%                   imagesc(pheromoneMat(1:50,1:50));
    end
    % reduce pheromone
    for u = 1:size(pheromoneMat,1)
        for r = 1:size(pheromoneMat,2)
            pheromoneMat(u,r)= pheromoneDecrement*pheromoneMat(u,r);
            if(pheromoneMat(u,r) < phimin)
                pheromoneMat(u,r) = phimin;
            end
        end
    end
end
 
%% end of iterations - ask all ants to drop and update map
for ant = 1:antCount
    if(antMouths(ant,1) ~= 0)
        foodMat(antLoc(ant,1),antLoc(ant,2)) = antMouths(ant,1);
    end
end
pause(.01);
cla;
plot(foodLoc(:,1),foodLoc(:,2),'go');
plot(antLoc(:,1),antLoc(:,2),'r*');

