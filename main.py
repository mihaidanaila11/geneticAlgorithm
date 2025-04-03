import math
import numpy as np
import bisect
from gui import windowInterface

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
        self.maxValue = float("-inf\n")

        def initPopulation(chromLen, chromNumber):
            bits = [str(int(x > 0.5)) for x in np.random.rand(chromNumber*chromLen)]

            population = []

            for i in range(chromNumber):
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
    
    def _evolveHelper(self, logFile, printInfo = False):
        # Afisam populatia initiala
        if printInfo:
            logFile.write("Populatia initiala\n")
            for i in range(self.populationSize):
                logFile.write(f"{i+1}: {self.population[i]} x={self.decodificare(self.population[i])} f={self.fitness(self.decodificare(self.population[i]))}\n")
            logFile.write("\n")

        # Lista de probabilitati
        selectionProbabilities = [self.getSelectionProbability(chromosome) for chromosome in self.population]
        if printInfo:
            for i in range(len(selectionProbabilities)):
                logFile.write(f"cromozom {i+1} probabilitate {selectionProbabilities[i]}\n")
        
        # Calculam intervalele pentru selectie
        cumulativeProb = [0]

        # Gasim cel mai bun cromozom pentru a-l baga mereu in selectie
        bestChromosomeIndex = 0

        for p in range(len(selectionProbabilities)):
            cumulativeProb.append(cumulativeProb[-1] + selectionProbabilities[p])

            if selectionProbabilities[p] > selectionProbabilities[bestChromosomeIndex]:
                bestChromosomeIndex = p

        if printInfo:
            logFile.write(f"Intervale probabilitati selectie\n{[float(x) for x in cumulativeProb]}\n") # Nu stiu cat de eficienta e conversia
            logFile.write(f"Cel mai bun cromozom este {bestChromosomeIndex+1}\n")       

        # Procesul de selectie
        randomPicks = np.random.rand(self.populationSize-1)

        selectii = [self.population[bestChromosomeIndex]]
        for u in randomPicks:
            # Gasesc din ce interal face parte u
            index = bisect.bisect_left(cumulativeProb, u) - 1

            # Il selectez
            if printInfo:
                logFile.write(f"u = {u} selectam cromozomul {index}\n")
            selectii.append(self.population[index])
        if printInfo:
            logFile.write("\n")

            logFile.write(f"Dupa selectie:\n")
            for i in range(self.populationSize):
                logFile.write(f"{i+1}: {selectii[i]} x={self.decodificare(selectii[i])} f={self.fitness(self.decodificare(selectii[i]))}\n")
            logFile.write("\n")

        # Procesul de incrucisare
        def crossover(chromosome1, chromosome2, index):
            if printInfo:
                logFile.write(f"{chromosome1} {chromosome2} punct {index}\n")

            result = (chromosome1[0:index] + chromosome2[index:], chromosome2[0:index] + chromosome1[index:])

            if printInfo:
                logFile.write(f"rezultat {result[0]} {result[1]}\n")

            return result
        
        randomPicks = np.random.rand(self.populationSize)
        puncteCrossover = np.random.randint(0, self.bitNumber, size=self.populationSize)

        indexChromosome1, indexChromosome2 = (None, None)

        if printInfo:
            logFile.write(f"Probabilitatea de incrucisare {self.crossoverProb}\n")

        for i in range(self.populationSize):
            if(indexChromosome1 and indexChromosome2):
                if printInfo:
                    logFile.write(f"recombinare dintre cromozomul {indexChromosome1} si cromozomul {indexChromosome2}:\n")

                selectii[indexChromosome1], selectii[indexChromosome2] = crossover(selectii[indexChromosome1], 
                                                                                   selectii[indexChromosome2],
                                                                                   puncteCrossover[indexChromosome1])
                indexChromosome1, indexChromosome2 = (None, None)


            if(randomPicks[i] < self.crossoverProb):
                if printInfo:
                    logFile.write(f"{i+1}: {selectii[i]} u={randomPicks[i]} < {self.crossoverProb} participa\n")
                
                if not indexChromosome1:
                    indexChromosome1 = i
                else:
                    indexChromosome2 = i
            else:
                if printInfo:
                    logFile.write(f"{i+1}: {selectii[i]} u={randomPicks[i]}\n")

        if printInfo:
            logFile.write("\n")

            logFile.write(f"Dupa recombinare:\n")
            for i in range(self.populationSize):
                logFile.write(f"{i+1}: {selectii[i]} x={self.decodificare(selectii[i])} f={self.fitness(self.decodificare(selectii[i]))}\n")
            logFile.write("\n")

        #mutatie
            logFile.write(f"Probabilitatea de mutatie penntru fiecare gena {self.mutationProb}\n")

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
                    logFile.write(f"A fost modificat cromozomul {chromosomeIndex+1}\n")
                selectii[chromosomeIndex] = "".join(modifiedChromosome)

        if printInfo:
            logFile.write(f"Dupa mutatie:\n")

        self.fitnessSum = 0
        for i in range(self.populationSize):
            fitness = self.fitness(self.decodificare(selectii[i]))
            self.fitnessSum += fitness
            
            # Pastram valoarea maxima a functiei gasita pana la generatia curenta
            if self.maxValue < fitness:
                self.maxValue = fitness

            if printInfo: 
                logFile.write(f"{i+1}: {selectii[i]} x={self.decodificare(selectii[i])} f={fitness}\n")
        if printInfo:       
            logFile.write("\n")
        
        self.population = selectii

        
    def evolve(self):
        evolutionLog = open("evolutie.txt", "w")
        self._evolveHelper(evolutionLog, printInfo = True)
        
        evolutionLog.write("Evolutia maximului:\n")
        evolutionLog.write(f"Max fitness: {self.maxValue}\n")
        evolutionLog.write(f"Fitness Mean: {self.fitnessSum / self.populationSize}\n\n")

        for _ in range(self.generationsNumber - 1):
            self._evolveHelper(evolutionLog, False)
            evolutionLog.write(f"Max fitness: {self.maxValue}\n")
            evolutionLog.write(f"Fitness Mean: {self.fitnessSum / self.populationSize}\n\n")

        evolutionLog.close()


# inputFile = open("input.txt", 'r')
# populationSize =        int(inputFile.readline().strip())
# left, right =           [float(x) for x in inputFile.readline().strip().split()]
# a,b,c =                 [int(x) for x in inputFile.readline().strip().split()]
# p =                     int(inputFile.readline().strip())
# crossoverProbability =  int(inputFile.readline().strip())
# mutationProbability =   int(inputFile.readline().strip())
# generationsNumber =     int(inputFile.readline().strip())
# inputFile.close()


gui = windowInterface(maxFunctie)
