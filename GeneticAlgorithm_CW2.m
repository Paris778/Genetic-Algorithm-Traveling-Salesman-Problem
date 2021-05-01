clc;
close all;
clear;

%% Load Dataset

load('xy.mat', 'xy');
dataset = xy;
disp(dataset);


%% Initialise Variables

%Iterations
iter = 2000;  
%Population Sizes
population_size = 100;

%Just to graph the performance of the algorithm per iteration
best_score = zeros(iter,1);
worst_score = zeros(iter,1);

[num_of_cities,dimentions] = size(dataset);

%Random sequence of cities
city_number = randperm(num_of_cities,100);

%Used for initial chromosome formation
temp_city_number = city_number;

disp("=======");
disp(dataset(98,:)); % Row , Column
disp(dataset(98,1));  %x

%% Calculate the distance between all possible city pairs
city_pair_distances = zeros(num_of_cities,num_of_cities);

for i=1 : num_of_cities
    for j=1 : num_of_cities
        city_pair_distances(i,j) = calculate_Distance(i,j,dataset);
    end 
end


%% Initialise Random Starting Chromosomes

population = zeros(population_size,num_of_cities);
%Make initial generation chromosomes
for i = 1:population_size
    %Random Chromosome
    temp_chromosome = randperm(num_of_cities,100); %%Permutation Encoding
    population(i,:) = temp_chromosome;
end


%% Always have an extra column at end for fitness calculation
population = [population zeros(population_size,1)];

%% Iterate over the Genetic Algorithm Process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
for k = 1:iter %How many times it will iterate
    
    % Replace distances with zeros
    population(:,num_of_cities+1) = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculate Fitness of each chromosome
    
    %For each member of the population, evauate their chromosome
     for i = 1:population_size
        %Loop over their chromosome
        for index = 1 : num_of_cities-1 %one less pair
            
            %distance = calculate_Distance(dataset, population, i,index); 
            distance = city_pair_distances(population(i,index),population(i,index+1));
            %%Add this result to the distance count
            population(i,num_of_cities+1) = distance + population(i,num_of_cities+1);
        end
        %Return back to starting point
        distance = city_pair_distances(population(i,num_of_cities),population(i,1));
        population(i,num_of_cities+1) = distance + population(i,num_of_cities+1);
        population(i,num_of_cities+1) = 1/population(i,num_of_cities+1);
     end

     best = max(population(:,num_of_cities+1));
     worst = min(population(:,num_of_cities+1));
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Elitist Selection, keep two best agents
    population = sortrows(population,num_of_cities+1);
    population_new = zeros(population_size,num_of_cities);
    
    %
    population_new(1:2,:) = population(population_size-1:population_size,1:num_of_cities);
    %population_new(1:2,:) = population(1:2,1:num_of_cities);
    population_new_num = 2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Make new population 
    while (population_new_num < population_size)
        %% Choose the two best chromosomes 
        weights= population(:,101)/sum(population(:,101));
        choice1 = Selection(weights);
        choice2 = Selection(weights);
        temp_chromosome_1 = population(choice1, 1:100);
        temp_chromosome_2 = population(choice2, 1:100);
        
        %% Crossover  between elite chromosomes
        if (rand < 0.75) %2-point crossover
              %[temp_chromosome_1,temp_chromosome_2] = crossover(temp_chromosome_1,temp_chromosome_2,num_of_cities);
        end
        
        %% Mutation prob 0.25
        
        %Chromosome 1 Mutation
        if (rand < 0.90)
             temp_chromosome_1 = mutate(temp_chromosome_1,num_of_cities);
        end
        
        %Chromosome 2 Mutation
        if (rand < 0.90)
             temp_chromosome_2 = mutate(temp_chromosome_2,num_of_cities);
        end
        
        %% Make a new population out of the mutated/crossed_Over chromosomes

            if (population_new_num < population_size)
                population_new_num = population_new_num + 1;
                population_new(population_new_num,:) = temp_chromosome_1;
            end

            if (population_new_num < population_size)
                population_new_num = population_new_num + 1;
                population_new(population_new_num,:) = temp_chromosome_2;
            end

    end
    
    
    disp("The champion of this iteration: ");
    disp(best); disp(worst);
    % Assign scores
    best_score(k,1) = 1/best;
    worst_score(k,1) = 1/worst;
    % Assign new population
    population(:,1:num_of_cities) = population_new;
    
end
toc
%% Find the best result
[best,best_index] = max(population(:,num_of_cities+1));
best = 1/best;
disp("===============================");
disp("Shortest Distance found: ");
disp(best);
disp("===============================");


%% Visualise 


optRoute = population(best_index,1:num_of_cities);
minDist = best;

 figure('Name','TSP_GA | Results','Numbertitle','off');
 subplot(2,2,1);
 pclr = ~get(0,'DefaultAxesColor');
 plot(xy(:,1),xy(:,2),'.','Color',pclr);
 title('City Locations');
 subplot(2,2,2);
 rte = optRoute([1:100 1]);
 plot(xy(rte,1),xy(rte,2),'r.-');
 title(sprintf('Total Distance = %1.4f',minDist));
 
 
 %% Visualise Performance %%%

% Top plot
figure(2);
ax1 = nexttile;
plot(ax1,1:iter,best_score,'g');
title(ax1,'Best Score Plot');
ylabel(ax1,'Distance');

% Bottom plot
ax2 = nexttile;
plot(ax2,1:iter,worst_score,'r');
title(ax2,'Worst Score Plot');
ylabel(ax2,'Distance');
 
 
 %% Calculte Fitness Function
 function distance = calculate_Distance(i,j, dataset)
 
        %Assign values to coordinates    
        
        x1 = dataset(i,1); % x1 value
        y1 = dataset(i,2); % y1 value
        %%%
        x2 = dataset(j,1); % x2 value
        y2 = dataset(j,2); % y2 value

        if(i == j)
            distance = 999999;
        else
            vector = [x1,y1;x2,y2];
            distance = pdist(vector, 'euclidean');
        end
 end
 
 
 
 %% Mutate Function 
function mutated_chromosome = mutate(temp_chromosome_1,num_of_cities)

              randomPoint = randi([1,num_of_cities],1);
              randomPoint_2 = randi([1,num_of_cities],1);

               %%%%%%%%%%%%%%%%
               %%Swap mutation                      
                  %temp_chromosome_1([randomPoint randomPoint_2]) = temp_chromosome_1([randomPoint_2 randomPoint]);

               %%%%%%%%%%%%%%%%
               %Slide Mutation                    
                  %temp_chromosome_1(randomPoint:randomPoint_2) = temp_chromosome_1([randomPoint+1:randomPoint_2 randomPoint]);

               %%%%%%%%%%%%%%%
               %Flip mutation
                  temp_chromosome_1(randomPoint:randomPoint_2) = temp_chromosome_1(randomPoint_2:-1:randomPoint);
                   
      mutated_chromosome = temp_chromosome_1;
end


%% Crossover Function
function [temp_chromosome_1_NEW,temp_chromosome_2_NEW] = crossover(temp_chromosome_1,temp_chromosome_2,num_of_cities)

    randomPoint = randi([2,num_of_cities-2],1);
    randomPoint_2 = randi([2,num_of_cities-2],1);

    %%Cross over from chrmsm_1 to chrmsm_2
        piece_1 = temp_chromosome_1(randomPoint:randomPoint_2);
        temp_chromosome_2_temp = setdiff(temp_chromosome_2,piece_1);
        %%Final:
        temp_chromosome_2_NEW = [temp_chromosome_2_temp(1:randomPoint) piece_1 temp_chromosome_2_temp(randomPoint+1:end)];
    
        %%Cross over from chrmsm_2 to chrmsm_1
        
        randomPoint = randi([2,num_of_cities-3],1);
        randomPoint_2 = randi([2,num_of_cities-3],1);
        
        piece_2 = temp_chromosome_2(randomPoint:randomPoint_2);
        temp_chromosome_1_temp = setdiff(temp_chromosome_1,piece_2);
        %%Final:
        temp_chromosome_1_NEW = [temp_chromosome_1_temp(1:randomPoint) piece_2 temp_chromosome_1_temp(randomPoint+1:end)];

end