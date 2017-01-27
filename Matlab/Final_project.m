clear all;

%go through and generate pseudo genes, at the set length
%Need even number of genes
geneCount =  50;
geneLength = 10;
for i = 1:geneCount
    genes{i} = '';
    for j = 1:geneLength
        choose = randi(4);
        if choose == 1
            genes{i} = strcat(genes{i},'A');
        end
        if choose == 2
            genes{i} = strcat(genes{i},'T');
        end
        if choose == 3
            genes{i} = strcat(genes{i},'C');
        end
        if choose == 4
            genes{i} = strcat(genes{i},'G');
        end
    end
end
%each position of these genes represents a 'state'
%Each state is an environmental factor (heat, cold, food...)
%generates a random weighting of each of these genes
for i = 1:geneCount
    geneWeight(i) = (randi(6) + 2) / 10;
end

%Each gene has a 'counter' state, such that it negtively affects the
%chances of surivial
%This generates an array, where the place in it is the state, and the
%number in it is its counter state

for i = 2:2:geneCount  
    counterLocation(i) = i - 1;
    counterLocation(i - 1) = i;
end

%chooses current environment states
maxNum = geneCount;
stateCount = randi(geneCount);
for i = 1:stateCount 
    if i == 1 
        states(i) = randi(maxNum);
    else
        found = false;
        state = randi(geneCount);
        for j = 1:size(states)
            if state == states(j)
                found = true;
            end
            if counterLocation(state) == states(j)
                found = true;
            end
        end
        if found ~= true
            states(i) = state;
        else
            states(i) = -1;
        end  
    end     
end

%Notes for implementation
%Create a cell array of species
%each species is a struct of the different members in that species
%loop through evolutions
%each time, it will check which genes the animal has, then, gene by gene,
%calculates first the effect of the gene on the animals survival, and then
%adjusts it for the genes weighting. Then it combines all of these together
%to get a survival number for this animal. This survival number affects the
%chances of the animal getting more resources. Each animal needs a certain
%number of resources or it dies, but it can take more than it needs. It
%then calculates how many resources each animal gets, and based on this,
%determines whether or not the animal dies. The animals that survive then
%have a chance of having offspring, and the more resource points they have,
%the more likely they are to have more offspring. These offspring then consitute
%the next generation.

%creates animals genome, put it into animal structs, and saves the indexes
%of the genes to these same structs then each of those into species

speciesCount = 20;
animalCount = 10;
genomeLength = 10;
modifier = 50;
for i = 1:speciesCount 
    animal = '';
    for j = 1:genomeLength
        num = randi(geneCount);
        geneLocs(j) = num;
        gene = genes{num};
        animal = strcat(animal, gene);     
    end
    for j = 1:animalCount
        animals(j).genome = animal;
        animals(j).genomeLength = genomeLength;
        animals(j).alive = true;
        animals(j).resourcePoints = 0;
        animals(j).currentPoints = 0;
        animals(j).modifier = modifier;
        animals(j).childNumber = 0;
        for k = 1:genomeLength
            animals(j).geneIndex(k) = geneLocs(k);
            animals(j).counterGenes(k) = counterLocation(geneLocs(k));
            animals(j).geneWeights(k) = geneWeight(geneLocs(k));    
        end
    end
    species{i} = animals;
end


generationNum = 100;
baseProb = 1;
pointsNeeded = 2;
generations{1} = species;
savedGenerations{1} = species;
reachedEnd = false;
for i = 1:generationNum
    generations{i} = species;
end
%Main evolution loop. 
for i = 1:generationNum
    assignedTempG = false;
    resourcePoints = 100000;
    changer = 0;
    if reachedEnd ~= true
        savedGenerations{i} = generations{i};
        speciesPoints = 0;
        %This loop calculates the raw points that every animal has by 
        %checking to see if the genes are affected by any of the environmental 
        %states, and adjusting the number appropriately. It also weights each
        %gene, and then adds this number together to produce the animals
        %score. All of the scores are then added together to get a total score
        %for all animals within a species, and once again to get them for all
        %of the species. 
            for j = 1:length(generations{i})
                animalPoints(j) = 0;
                for k = 1:length(generations{i}{j})
                    totalPoints = 0;
                        totalWeight = 0;
                        for l = 1:generations{i}{j}(k).genomeLength
                            totalWeight = totalWeight + generations{i}{j}(k).geneWeights(l);
                        end
                        for l = 1:generations{i}{j}(k).genomeLength
                            index = generations{i}{j}(k).geneIndex(l);
                            counterIndex = generations{i}{j}(k).counterGenes(l);
                            weight = generations{i}{j}(k).geneWeights(l) / totalWeight;
                            found1 = false;
                            found2 = false;
                            for m = 1:stateCount
                                if index == states(m)
                                    found1 = true;
                                end
                                if counterIndex == states(m)
                                    found2 = true;
                                end
                            end
                            if found1 == true
                                change = baseProb + generations{i}{j}(k).modifier;
                                points = weight * change;
                                totalPoints = totalPoints + points;
                            elseif found2 == true
                                change = baseProb - generations{i}{j}(k).modifier;
                                points = weight * change;
                                totalPoints = totalPoints + points;
                            else 
                                change = baseProb;
                                points = weight * change;
                                totalPoints = totalPoints + points;
                            end
                        end
                    generations{i}{j}(k).currentPoints = totalPoints;
                    animalPoints(j) =  animalPoints(j) + totalPoints;
                end
                    speciesPoints = speciesPoints + animalPoints(j);
            end
            clear speciesIsDead;
            for j = 1:length(generations{i})
                speciesIsDead(j) = false;
            end
            %This weights the total number of points to get the number of resource
            %points that each animal will recieve. It then checks to see if the
            %animals have the minimum number of points needed to survive, and if
            %they don't they die. It also checks to see if an entire species
            %died, and if it did, it removes it from the current species
            pointsUsed = 0;
            for j = 1:length(generations{i})
                isAlive = false;
                for k = 1:length(generations{i}{j})
                    percentPoints = generations{i}{j}(k).currentPoints / speciesPoints;
                    seedNum = floor(percentPoints * 10000);
                    if seedNum > 0
                        pointsForAnimal = randi(seedNum);
                    else 
                        pointsForAnimal = 0;
                    end
                    resourcePoints = resourcePoints - pointsForAnimal;
                    if resourcePoints > 0
                        generations{i}{j}(k).resourcePoints = pointsForAnimal;
                    else 
                        generations{i}{j}(k).resourcePoints = 0;
                    end
                    if generations{i}{j}(k).resourcePoints < pointsNeeded
                        generations{i}{j}(k).alive = false;
                    elseif isAlive == false
                        isAlive = true;
                    end           
                end
                if isAlive == false 
                    speciesIsDead(j) = true;
                end
            end
            aliveCounter = 0;
            for j = 1:length(speciesIsDead)
                if speciesIsDead(j) == false
                    aliveCounter = aliveCounter + 1;
                    speciesIndex(j) = j;
                else 
                    speciesIndex(j) = -1;
                end
            end
            myCounter = 1;
            clear currentSpecies;
            for j = 1:length(generations{i})
                if speciesIndex(j) ~= -1
                    currentSpecies{myCounter} = generations{i}{speciesIndex(j)};
                    myCounter = myCounter + 1;
                end
            end
            %creates a new generation of animals, getting rid of any
            %species that died off. 
            birtherLocation(1) = 0;      
            if myCounter ~= 1
                for j = 1:length(currentSpecies)
                    hasChildren = false;
                    mySize = 1;
                    birthers = 0;
                    for k = 1:length(currentSpecies{j})
                        reproductionNumber = currentSpecies{j}(k).resourcePoints - (pointsNeeded);
                        if currentSpecies{j}(k).alive == true
                           childNumber =   floor((rand) * 5 * sqrt(sqrt(reproductionNumber)));
                           currentSpecies{j}(k).childNumber = childNumber;
                           if childNumber > 0
                               birthers = birthers + 1;
                               birtherLocation(birthers) = k;
                               hasChildren = true; 
                           end           
                        end     
                    end
                    counter = 1;    
                    for k = 1:birthers               
                            for l = 1:currentSpecies{j}(birtherLocation(k)).childNumber
                                newAnimals(counter).genome = currentSpecies{j}(birtherLocation(k)).genome;
                                newAnimals(counter).genomeLength = genomeLength;
                                newAnimals(counter).alive = true;
                                newAnimals(counter).resourcePoints = 0;
                                newAnimals(counter).currentPoints = 0;
                                newAnimals(counter).modifier = modifier;
                                newAnimals(counter).childNumber = 0;
                                for m = 1:genomeLength
                                    newAnimals(counter).geneIndex(m) = currentSpecies{j}(birtherLocation(k)).geneIndex(m);
                                    newAnimals(counter).counterGenes(m) = counterLocation(currentSpecies{j}(birtherLocation(k)).geneIndex(m));
                                    newAnimals(counter).geneWeights(m) = geneWeight(currentSpecies{j}(birtherLocation(k)).geneIndex(m));
                                end
                                counter = counter + 1;  
                            end
                    end       
                    clear birtherLocation;
                    clear birthers;
                    if hasChildren == true
                        tempGeneration{j - changer} = newAnimals;
                        assignedTempG = true;
                        clear newAnimals;
                    else 
                        changer = changer + 1;
                        tempGenerations{j} = [];
                    end     
                    clear speciesIndex;
                end    
            else 
                reachedEnd = true;
            end
            if reachedEnd ~= true && assignedTempG == true
                clear generations{i + 1}; 
                generations{i + 1} = tempGeneration;
                clear tempGeneration
            elseif assignedTempG == false
                reachedEnd = true;
                clear generations{i + 1};
                clear tempGeneration;
            else 
                clear generations{i + 1};
                clear tempGeneration;
            end
    else
        generations{i} = [];
    end
 
end
%reorganizes data structures so that it can be graphed more easily
time = 1:generationNum;
max = -1;
for i = 1:generationNum 
    if length(generations{i}) > max
        max = length(generations{i});
    end
end
for i = 1:length(generations{1})
    trackGenes{i} = generations{1}{i}(1).genome;
end
    
clear animals;
for i = 1:generationNum 
    for j = 1:length(trackGenes)
        found = false;
        for k = 1:length(generations{i})
            if trackGenes{j} == generations{i}{k}(1).genome
                found = true;
                index = k;
            end
        end
        if found == true
            animals{j}(i).size = length(generations{i}{index});
        else 
            animals{j}(i).size = 0;
        end
    end
end
for i = 1:length(animals)
    for j = 1:length(animals{i})
        plotter(j) = animals{i}(j).size;
    end
    hold on;
    plot(time, plotter, 'Color',rand(1,3));
    xlabel('Generation');
    ylabel('Population of species');
    title('Simulation of the laws of survival of the fittest');
end

























    
    
   














