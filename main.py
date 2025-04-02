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
        self.mutationProb = mutationProb

        # Numarul de generatii
        self.generationsNumber = generationsNumber

        # Populatia

        # Numarul de biti necesari pentru codificare
        self.bitNumber = math.ceil(abs(math.log2((right-left)*(pow(10,p)))))

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
    
    def evolve(self):
        
        

        # Lista de probabilitati
        selectionProbabilities = [self.getSelectionProbability(chromosome) for chromosome in self.population]
        
        # Calculam intervalele pentru selectie
        cumulativeProb = [0]
        for p in selectionProbabilities:
            cumulativeProb.append(cumulativeProb[-1] + p)


        # Procesul de selectie
        randomPicks = np.random.rand(self.populationSize)

        selectii = []
        for u in randomPicks:
            # Gasesc din ce interal face parte u
            index = bisect.bisect_left(cumulativeProb, u) - 1

            # Il selectez
            print(f"u = {u} selectam cromozomul {index}")
            selectii.append(self.population[index])

        print(f"Dupa selectie:")
        for i in range(self.populationSize):
            print(f"{i+1}: {selectii[i]} x={self.decodificare(selectii[i])} f={self.fitness(self.decodificare(selectii[i]))}")

        # Procesul de incrucisare
        def crossover(chromosome1, chromosome2, index):
            print(f"{chromosome1} {chromosome2} punct {index}")

            result = (chromosome1[0:index] + chromosome2[index:], chromosome2[0:index] + chromosome1[index:])

            print(f"rezultat {result[0]} {result[1]}")

            return result
        
        randomPicks = np.random.rand(self.populationSize)
        puncteCrossover = np.random.randint(0, self.bitNumber, size=self.population)

        indexChromosome1, indexChromosome2 = (None, None)

        print(f"Probabilitatea de incrucisare {self.crossoverProb}")
        for i in range(self.populationSize):
            if(indexChromosome1 and indexChromosome2):
                print(f"recombinare dintre cromozomul {indexChromosome1} si cromozomul {indexChromosome2}:")
                selectii[indexChromosome1], selectii[indexChromosome2] = crossover(selectii[indexChromosome1], 
                                                                                   selectii[indexChromosome2],
                                                                                   puncteCrossover[indexChromosome1])
                indexChromosome1, indexChromosome2 = (None, None)


            if(randomPicks[i] < self.crossoverProb):
                print(f"{i+1}: {selectii[i]} u={randomPicks[i]} < {self.crossoverProb}")
                
                if not indexChromosome1:
                    indexChromosome1 = i
                else:
                    indexChromosome2 = i
            else:
                print(f"{i+1}: {selectii[i]} u={randomPicks[i]}")

        #mutatie
        #adaugam elementele elitiste

        



    


algGen = maxFunctie(20, -1, 2, -1, 1, 2, 6, 25, 1, 50)

populatie = [int(binary, 2) for binary in algGen.population]

algGen.evolve()
