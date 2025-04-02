import math
import numpy as np
import bisect

class maxFunctie:
    def decToBin(self, x):
        ans = ""
        while x:
            ans = str(x & 1) + ans
            x = x >> 1

        return ans  


    def fitness(self, x):
        return self.a * x * x + self.b * x + self.c
    def __init__(self, populationSize, left, right, a, b, c, p, crossoverProb, mutationProb,
                 generationsNumber):
        # Dimensiunea populatiei
        self.populationSize = populationSize

        # Domeniul de definite al functiei
        self.left = left
        self.right = right

        # Parametrii functiei
        self.a = a
        self.b = b
        self.c = c

        # Precizia cu care se lucreaza
        self.p = p

        # Probabilitatea de incrucisare
        self.crossoverProb = crossoverProb / 100

        # Probabilitatea de mutatie
        self.mutationProb = mutationProb / 100

        # Numarul de generatii
        self.generationsNumber = generationsNumber

        # Populatia

        # Numarul de biti necesari pentru codificare
        self.bitNumber = math.ceil(abs(math.log2((right-left)*(pow(10,p)))))

        # Maximul functiei
        self.maxValue = float("-inf")

        def initPopulation(chromLen, chronNumber):
            bits = [str(int(x > 0.5)) for x in np.random.rand(chronNumber*chromLen)]

            population = []

            for i in range(chronNumber):
                population.append("".join(bits[i * chromLen : i * chromLen + chromLen]))

            return population
        
        self.population = initPopulation(self.bitNumber, populationSize)

        self.disc = (right-left)/pow(2, self.bitNumber)

    def codificare(self, x):
        k = self.decToBin(int((x- self.left)//self.disc))

        k = "0" * (self.bitNumber - len(k)) + k

        return k

    def decodificare(self, x):
        x = int(x, 2)

        k = self.left + x * self.disc

        return k
    
    def getTotalPerformance(self):
        return np.sum([self.fitness(self.decodificare(chromosome)) for chromosome in self.population])
    
    def getSelectionProbability(self, chromosome):
        return self.fitness(self.decodificare(chromosome)) / self.getTotalPerformance()
    
    def _evolveHelper(self, printInfo = False):
        # Afisam populatia initiala
        if printInfo:
            print("Populatia initiala")
            for i in range(self.populationSize):
                print(f"{i+1}: {self.population[i]} x={self.decodificare(self.population[i])} f={self.fitness(self.decodificare(self.population[i]))}")
            print()

        # Lista de probabilitati
        selectionProbabilities = [self.getSelectionProbability(chromosome) for chromosome in self.population]
        
        # Calculam intervalele pentru selectie
        cumulativeProb = [0]

        # Gasim cel mai bun cromozom pentru a-l baga mereu in selectie
        bestChromosomeIndex = 0

        for p in range(len(selectionProbabilities)):
            cumulativeProb.append(cumulativeProb[-1] + selectionProbabilities[p])

            if selectionProbabilities[p] > selectionProbabilities[bestChromosomeIndex]:
                bestChromosomeIndex = p

        if printInfo:
            print(f"Cel mai bun cromozom este {bestChromosomeIndex+1}")       

        # Procesul de selectie
        randomPicks = np.random.rand(self.populationSize-1)

        selectii = [self.population[bestChromosomeIndex]]
        for u in randomPicks:
            # Gasesc din ce interal face parte u
            index = bisect.bisect_left(cumulativeProb, u) - 1

            # Il selectez
            if printInfo:
                print(f"u = {u} selectam cromozomul {index}")
            selectii.append(self.population[index])
        if printInfo:
            print()

            print(f"Dupa selectie:")
            for i in range(self.populationSize):
                print(f"{i+1}: {selectii[i]} x={self.decodificare(selectii[i])} f={self.fitness(self.decodificare(selectii[i]))}")
            print()

        # Procesul de incrucisare
        def crossover(chromosome1, chromosome2, index):
            if printInfo:
                print(f"{chromosome1} {chromosome2} punct {index}")

            result = (chromosome1[0:index] + chromosome2[index:], chromosome2[0:index] + chromosome1[index:])

            if printInfo:
                print(f"rezultat {result[0]} {result[1]}")

            return result
        
        randomPicks = np.random.rand(self.populationSize)
        puncteCrossover = np.random.randint(0, self.bitNumber, size=self.populationSize)

        indexChromosome1, indexChromosome2 = (None, None)

        if printInfo:
            print(f"Probabilitatea de incrucisare {self.crossoverProb}")

        for i in range(self.populationSize):
            if(indexChromosome1 and indexChromosome2):
                if printInfo:
                    print(f"recombinare dintre cromozomul {indexChromosome1} si cromozomul {indexChromosome2}:")

                selectii[indexChromosome1], selectii[indexChromosome2] = crossover(selectii[indexChromosome1], 
                                                                                   selectii[indexChromosome2],
                                                                                   puncteCrossover[indexChromosome1])
                indexChromosome1, indexChromosome2 = (None, None)


            if(randomPicks[i] < self.crossoverProb):
                if printInfo:
                    print(f"{i+1}: {selectii[i]} u={randomPicks[i]} < {self.crossoverProb}")
                
                if not indexChromosome1:
                    indexChromosome1 = i
                else:
                    indexChromosome2 = i
            else:
                if printInfo:
                    print(f"{i+1}: {selectii[i]} u={randomPicks[i]}")

        if printInfo:
            print()

            print(f"Dupa recombinare:")
            for i in range(self.populationSize):
                print(f"{i+1}: {selectii[i]} x={self.decodificare(selectii[i])} f={self.fitness(self.decodificare(selectii[i]))}")
            print()

        #mutatie
            print(f"Probabilitatea de mutatie penntru fiecare gena {self.mutationProb}")

        # Varianta 2 de recombinare
        randomPicks = np.random.rand(self.populationSize * self.bitNumber)
        for chromosomeIndex in range(len(selectii)):
            chromosome = selectii[chromosomeIndex]
            modifiedChromosome = list(chromosome)
            modifiedFlag = False
            for i in range (len(chromosome)):
                u = randomPicks[chromosomeIndex * self.bitNumber + i]
                if u < self.mutationProb:
                    modifiedChromosome[i] = '0' if modifiedChromosome[i] == '1' else '1'
                    modifiedFlag = True

            if modifiedFlag:
                if printInfo:
                    print(f"A fost modificat cromozomul {chromosomeIndex+1}")
                selectii[chromosomeIndex] = "".join(modifiedChromosome)

        if printInfo:
            print(f"Dupa mutatie:")
        for i in range(self.populationSize):
            fitness = self.fitness(self.decodificare(selectii[i]))
            
            # Pastram valoarea maxima a functiei gasita pana la generatia curenta
            if self.maxValue < fitness:
                self.maxValue = fitness

            if printInfo: 
                print(f"{i+1}: {selectii[i]} x={self.decodificare(selectii[i])} f={fitness}")
        if printInfo:       
            print()
        
        self.population = selectii

        
    def evolve(self):
        self._evolveHelper(printInfo = True)
        
        print("Evolutia maximului:")
        print(self.maxValue)

        for _ in range(self.generationsNumber - 1):
            self._evolveHelper(printInfo = False)
            print(self.maxValue)



    

algGen = maxFunctie(20, -1, 2, -1, 1, 2, 6, 25, 1, 50)

populatie = [int(binary, 2) for binary in algGen.population]

algGen.evolve()
