import random as r
import os as o
import csv as c
from datetime import date as d


def clear():
    print("\033[H\033[2J", end="")


def continueText():
    input('Press "ENTER" to continue\n>>> ')


def createFolders():
    folders = ["Save", "Load", "Results"]
    for folder in folders:
        if not o.path.exists(folder):
            o.makedirs(folder)


class GeneticEntity:
    def __init__(self, symbol, descriptors):
        self._symbol = symbol
        self._descriptors = descriptors

    def getSymbol(self):
        return self._symbol

    def getDescriptors(self):
        return self._descriptors

    def _getGeneticSequence(self, geneticID):
        return "".join(self._symbol.upper() if int(digit) else self._symbol.lower() for digit in geneticID)

    def getGeneticSequence(self, geneticID):
        return self._getGeneticSequence(geneticID)


class Phenotype(GeneticEntity):
    def __init__(self, trait, symbol, descriptors):
        super().__init__(symbol, descriptors)
        self.__trait = trait

    def getTrait(self):
        return self.__trait

    def createGenotype(self, configuration):
        return Genotype(configuration, self.getSymbol(), self.getDescriptors())


class Genotype(GeneticEntity):
    def __init__(self, geneticID, symbol, descriptors):
        super().__init__(symbol, descriptors)
        self.__geneticID = geneticID

    def getGeneticID(self):
        return self.__geneticID

    def setGeneticID(self, value):
        self.__geneticID = value

    def getGeneticSequenceID(self):
        return self._getGeneticSequence(self.__geneticID)

    def getTrait(self):
        if self.__geneticID in ("11", "10"):
            traitValue = 1
        elif self.__geneticID == "00":
            traitValue = 0
        else:
            traitValue = -1

        if 0 <= traitValue < len(self._descriptors):
            return self._descriptors[traitValue]
        return "?"

    def display(self):
        print(f"{self.getGeneticSequenceID()} ({self.getTrait()})")


class Genome:
    def __init__(self, genotypes, phenotypes):
        self.__genotypes = genotypes
        self.__phenotypes = phenotypes

    def getGenotypes(self):
        return self.__genotypes

    def getPhenotypes(self):
        return self.__phenotypes

    def display(self):
        for genotype in self.__genotypes:
            genotype.display()

    @staticmethod
    def _normaliseID(geneticID):
        return "10" if geneticID == "01" else geneticID

    def crossover(self, parent):
        childGeneticIDs = [r.choice([self._normaliseID(a1 + a2)
                           for a1 in self.__genotypes[i].getGeneticID()
                           for a2 in parent.getGenotypes()[i].getGeneticID()])
                           for i in range(len(self.__genotypes))]

        childGenotypes = [self.__phenotypes[i].createGenotype(childGeneticIDs[i]) for i in range(len(self.__genotypes))]
        return Genome(childGenotypes, self.__phenotypes)


class Simulation:
    def __init__(self, mutationRate, phenotypes, elitismRate, initialPopulationSize, simName, genNum):
        self.__mutationRate = mutationRate
        self.__simName = simName
        self.__phenotypes = phenotypes
        self.__elitismRate = elitismRate
        self.__initialPopulationSize = initialPopulationSize
        self.__population = self.__initialisePopulation(initialPopulationSize)
        self.__generation = 0
        self.__genNum = genNum
        self.__data = []
        self.__parameters = {"mutationRate": mutationRate, "elitismRate": elitismRate, "initialPopulationSize": initialPopulationSize, "genNum": genNum}

    def getMutationRate(self):
        return self.__mutationRate

    def getSimName(self):
        return self.__simName

    def getPhenotypes(self):
        return self.__phenotypes

    def getElitismRate(self):
        return self.__elitismRate

    def getInitialPopulationSize(self):
        return self.__initialPopulationSize

    def getPopulation(self):
        return self.__population

    def getGeneration(self):
        return self.__generation

    def getGenNum(self):
        return self.__genNum

    def getData(self):
        return self.__data

    def getParameters(self):
        return self.__parameters

    def setPopulation(self, population):
        self.__population = population

    def __initialisePopulation(self, size):
        population = []
        for i in range(size):
            genotypes = [phenotype.createGenotype(self.__randomGeneticID()) for phenotype in self.__phenotypes]
            population.append(Genome(genotypes, self.__phenotypes))
        return population

    def __randomGeneticID(self):
        geneticID = "01"
        while geneticID == "01":
            geneticID = "".join(r.choice(["0", "1"]) for i in range(2))
        return geneticID

    def runSimulation(self, generations):
        for i in range(generations):
            self.__generation += 1
            newPopulation = self.__performGeneticAlgorithm()
            self.__updatePopulation(newPopulation)
            self.__collectData()
            self.__displayStatistics()
        self.__queryInit()
        self.__queryResults()

    def __performGeneticAlgorithm(self):
        newPopulation = []
        populationSize = len(self.__population)

        eliteCount = int(populationSize * self.__elitismRate)
        sortedPopulation = sorted(
            self.__population, key=self.__calculateFitness, reverse=True)
        newPopulation.extend(sortedPopulation[:eliteCount])

        fitnessSum = sum(self.__calculateFitness(g) for g in self.__population)
        weights = [self.__calculateFitness(g)/fitnessSum if fitnessSum != 0 else 1/populationSize for g in self.__population]

        for i in range(populationSize):
            parents = r.choices(self.__population, weights=weights, k=2)
            child = parents[0].crossover(parents[1])
            self.__mutate(child)
            newPopulation.append(child)

        return newPopulation

    def __mutate(self, genome):
        for genotype in genome.getGenotypes():
            newID = ""
            for allele in genotype.getGeneticID():
                if r.random() < self.__mutationRate:
                    newID += "1" if allele == "0" else "0"
                else:
                    newID += allele
            genotype.setGeneticID(Genome._normaliseID(newID))

    def __calculateFitness(self, genome):
        fitness = 0
        traitWeights = {"Speed": 1.0, "Metabolism": 1.2, "Bravado": 1.5}
        for i in range(len(genome.getGenotypes())):
            geneticID = genome.getGenotypes()[i].getGeneticID()
            fitness += int(geneticID, 2) * traitWeights.get(self.__phenotypes[i].getTrait(), 1.0)
        return fitness

    def __updatePopulation(self, newPopulation):
        self.__population = newPopulation
        survivalThreshold = 2
        self.__population = [genome for genome in self.__population if self.__calculateFitness(
            genome) > survivalThreshold]

    def __collectData(self):
        populationSize = len(self.__population)
        totalFitness = sum(self.__calculateFitness(genome) for genome in self.__population)
        averageFitness = totalFitness / populationSize if populationSize != 0 else 0

        genotypeDistribution = {}
        for genome in self.__population:
            for genotype in genome.getGenotypes():
                geneticID = genotype.getGeneticID()
                genotypeDistribution[geneticID] = genotypeDistribution.get(geneticID, 0) + 1

        phenotypeFrequency = {phenotype.getTrait(): {} for phenotype in self.__phenotypes}
        for genome in self.__population:
            for i in range(len(genome.getGenotypes())):
                genotype = genome.getGenotypes()[i]
                phenotype = self.__phenotypes[i]
                trait = genotype.getTrait()
                phenotypeFrequency[phenotype.getTrait()][trait] = phenotypeFrequency[phenotype.getTrait()].get(trait, 0) + 1

        self.__data.append({"generation": self.__generation, "populationSize": populationSize, "averageFitness": averageFitness, "genotypeDistribution": genotypeDistribution, "phenotypeFrequency": phenotypeFrequency})

    def __displayStatistics(self):
        clear()
        latestData = self.__data[-1]
        print(f"Generation {self.__generation}:")
        print(f"Population Size: {latestData['populationSize']}")
        print(f"Average Fitness: {latestData['averageFitness']:.2f}")
        print(f"Genotype Distribution: {latestData['genotypeDistribution']}")
        print("Phenotype Frequency:")
        for trait, frequencies in latestData['phenotypeFrequency'].items():
            print(f" {trait}: {frequencies}")
        continueText()

    def __queryInit(self):
        clear()
        ans = input("Save initial parameters? (Y/N): ")
        while ans.upper() not in ["Y", "N"]:
            print("Please enter Y or N.")
            ans = input("Save initial parameters? (Y/N): ")

        if ans.upper() == "Y":
            self.__saveInit()
        continueText()

    def __saveInit(self):
        with open(f"./Save/{self.__simName}.evsm", mode="w", newline="") as file:
            writer = c.writer(file)
            writer.writerow(["Parameter", "Value"])
            for key, value in self.__parameters.items():
                writer.writerow([key, value])
        print(f"Parameters saved to ./Save/{self.__simName}.evsm")

    def __queryResults(self):
        clear()
        ans = input("Save results? (Y/N): ")
        while ans.upper() not in ["Y", "N"]:
            print("Please enter Y or N.")
            ans = input("Save initial parameters? (Y/N): ")

        if ans.upper() == "Y":
            self.__saveResults()
        continueText()

    def __saveResults(self):
        filename = f"./Results/results-{self.__simName}-{d.today().strftime('%d-%m-%Y')}.txt"
        with open(filename, mode="w", newline="") as file:
            file.write("Simulation Results\n\n")
            for record in self.__data:
                file.write(f"Generation {record['generation']}:\n")
                file.write(f"Population Size: {record['populationSize']}\n")
                file.write(f"Average Fitness: {record['averageFitness']:.2f}\n")
                file.write(f"Genotype Distribution: {record['genotypeDistribution']}\n")
                for trait, frequencies in record['phenotypeFrequency'].items():
                    file.write(f" {trait}: {frequencies}\n")
                file.write("\n")
        print(f"Results saved to {filename}")


class SimulationMenu:
    def __init__(self):
        self.__running = True
        self.__simName = ""
        self.__mutationRate = 0.0
        self.__elitismRate = 0.0
        self.__initialPopulationSize = 0
        self.__genNum = 0

    def displayMenu(self):
        while self.__running:
            clear()
            print("EvoSim\n")
            print("1. Start Simulation")
            print("2. Load Simulation")
            print("3. Tutorial")
            print("4. Learn")
            print("5. Exit")
            choice = input(">>> ")
            self.__handleChoice(choice)

    def __handleChoice(self, choice):
        if choice == "1":
            self.__startSimulation()
        elif choice == "2":
            self.__loadSimulation()
        elif choice == "3":
            self.showTutorials()
        elif choice == "4":
            self.showLearn()
        elif choice == "5":
            self.exitProgram()
        else:
            print("Please choose a number between 1 and 5.")
            continueText()

    def __startSimulation(self):
        self.__setParameters()
        phenotypes = self.__initialisePhenotypes()
        simulation = Simulation(self.__mutationRate, phenotypes, self.__elitismRate, self.__initialPopulationSize, self.__simName, self.__genNum)
        self.__displayInitialPopulation(simulation)
        simulation.runSimulation(self.__genNum)

    def __setParameters(self):
        clear()
        print("Set Simulation Parameters:")
        self.__simName = input("Enter simulation name (e.g: example1): ")
        self.__mutationRate = self.__getFloatInput("Enter mutation rate (e.g: 0.01): ", 0.0, 1.0)
        self.__elitismRate = self.__getFloatInput("Enter elitism rate (e.g: 0.4): ", 0.0, 1.0)
        self.__initialPopulationSize = self.__getIntInput("Enter initial population size (e.g: 12): ", 2, 20)
        self.__genNum = self.__getIntInput("Enter number of generations (e.g: 7): ", 1, 14)
        continueText()

    def __getFloatInput(self, prompt, minValue, maxValue):
        while True:
            try:
                value = float(input(prompt))
                if minValue <= value <= maxValue:
                    return value
                print(f"Please enter a value between {minValue} and {maxValue}.")
            except ValueError:
                print("Invalid input. Please enter a number.")

    def __getIntInput(self, prompt, minValue, maxValue):
        while True:
            try:
                value = int(input(prompt))
                if minValue <= value <= maxValue:
                    return value
                print(f"Please enter an integer between {minValue} and {maxValue}.")
            except ValueError:
                print("Invalid input. Please enter an integer.")

    def __initialisePhenotypes(self):
        speed = Phenotype("Speed", "s", ["Slow", "Fast"])
        metabolism = Phenotype("Metabolism", "m", ["Slow Metabolism", "Fast Metabolism"])
        bravado = Phenotype("Bravado", "b", ["Timid", "Courageous"])
        return [speed, metabolism, bravado]

    def __displayInitialPopulation(self, simulation):
        clear()
        print("Initial Population:")
        for genome in simulation.getPopulation():
            genome.display()
            print()
        continueText()

    def __loadParametersFromCSV(self, filename):
        parameters = {}
        with open(filename, mode="r") as file:
            reader = c.reader(file)
            next(reader)
            for row in reader:
                parameters[row[0]] = row[1]
        return parameters

    def __loadSimulation(self):
        clear()
        print("Load Simulation")
        files = [f for f in o.listdir("./Load") if f.endswith(".evsm")]
        if not files:
            print("No .evsm files found in the Load folder.")
            continueText()

        print("Available .evsm files:")
        i = 1
        for file in files:
            print(f"{i}. {file}")
            i += 1

        choice = self.__getIntInput("Enter file number to load: ", 1, len(files))
        filename = files[choice - 1]
        filepath = f"./Load/{filename}"

        parameters = self.__loadParametersFromCSV(filepath)
        self.__mutationRate = float(parameters["mutationRate"])
        self.__elitismRate = float(parameters["elitismRate"])
        self.__initialPopulationSize = int(parameters["initialPopulationSize"])
        self.__genNum = int(parameters["genNum"])
        self.__simName = filename.split(".")[0]

        phenotypes = self.__initialisePhenotypes()
        simulation = Simulation(self.__mutationRate, phenotypes, self.__elitismRate, self.__initialPopulationSize, self.__simName, self.__genNum)
        self.__displayInitialPopulation(simulation)
        simulation.runSimulation(self.__genNum)

    def showTutorials(self):
        clear()
        print("Tutorial:")
        continueText()

    def showLearn(self):
        clear()
        print("Learn:")
        print("The predominant current-day meaning of genotype is some relevant part of the DNA passed to the organism by its parents. The phenotype is the physical and behavioral traits of the organism, for example, size and shape, metabolic activities, and patterns of movement. The distinction between them is especially important in evolutionary theory, where the survival and mating of organisms depends on their traits, but it is the DNA, held to be unaffected by the development of the traits over the life course, that is transmitted to the next generation.")
        continueText()

    def exitProgram(self):
        clear()
        print("Exiting...")
        self.__running = False


if __name__ == "__main__":
    createFolders()
    menu = SimulationMenu()
    menu.displayMenu()
