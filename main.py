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
    def __init__(self, symbol: chr, descriptors: list):
        self.symbol = symbol
        self.descriptors = descriptors

    def getGeneticSequence(self, geneticID: str):
        return "".join(self.symbol.upper() if int(digit) else self.symbol.lower() for digit in geneticID)


class Phenotype(GeneticEntity):
    def __init__(self, trait: str, symbol: chr, descriptors: list):
        super().__init__(symbol, descriptors)
        self.trait = trait

    def createGenotype(self, configuration: str):
        return Genotype(configuration, self.symbol, self.descriptors)


class Genotype(GeneticEntity):
    def __init__(self, geneticID: str, symbol: chr, descriptors: list):
        super().__init__(symbol, descriptors)
        self.geneticID = geneticID

    def getGeneticID(self):
        return self.getGeneticSequence(self.geneticID)

    def getTrait(self):
        if self.geneticID == "11" or self.geneticID == "10":
            traitValue = 1  # Stronger descriptor
        elif self.geneticID == "00":
            traitValue = 0  # Weaker descriptor
        else:
            traitValue = -1  # Undefined case

        # Ensure traitValue is valid (DEBUG)
        if 0 <= traitValue < len(self.descriptors):
            return self.descriptors[traitValue]
        else:
            print(self.descriptors)
            return "?"
    
    def display(self):
        print(f"{self.getGeneticID()} ({self.getTrait()})")


class Genome:
    def __init__(self, genotypes: list, phenotypes: list):
        self.genotypes = genotypes
        self.phenotypes = phenotypes

    def display(self):
        for genotype in self.genotypes:
            genotype.display()

    def normaliseID(self, geneticID: str):
        return "10" if geneticID == "01" else geneticID

    def crossover(self, parent):
        childGeneticIDs = [r.choice([self.normaliseID(allele1 + allele2) for allele1 in self.genotypes[element].geneticID for allele2 in parent.genotypes[element].geneticID]) for element in range(len(self.genotypes))]
        childGenotypes = [self.phenotypes[element].createGenotype(childGeneticIDs[element]) for element in range(len(self.genotypes))]
        return Genome(childGenotypes, self.phenotypes)

class Simulation:
    def __init__(self, mutationRate: float, phenotypes: list, elitismRate: float, initialPopulationSize: int, simName: str, genNum: int):
        self.mutationRate = mutationRate
        self.simName = simName 
        self.phenotypes = phenotypes
        self.elitismRate = elitismRate
        self.population = self.initialisePopulation(initialPopulationSize) 
        self.generation = 0
        self.genNum = genNum
        self.data = []
        self.parameters = {
            "Mutation Rate": mutationRate, 
            "Elitism Rate": elitismRate,
            "Initial Population Size": initialPopulationSize,
            "Number of Generations": genNum
        }
    
    def initialisePopulation(self, size: int): 
        population = [] 
        for i in range(size): 
            genotypes = [phenotype.createGenotype(self.randomGeneticID()) for phenotype in self.phenotypes]
            genome = Genome(genotypes, self.phenotypes) 
            population.append(genome) 
        return population 
    
    def randomGeneticID(self): 
        geneticID = "01"
        while geneticID == "01": 
            geneticID = "".join(r.choice(["0", "1"]) for i in range(2)) 
        return geneticID
    
    def runSimulation(self, generations):
        for i in range(generations):
            self.generation += 1
            newPopulation = self.performGeneticAlgorithm()
            self.updatePopulation(newPopulation)
            self.collectData()
            self.displayStatistics()
        self.queryInit()
        self.queryResults()

    def queryInit(self):
        clear()
        ans = input("Save initial parameters? (Y/N): ")
        while ans.upper() not in ["Y", "N"]:
            print("Please enter Y or N.")
            ans = input("Save initial parameters? (Y/N): ")
        
        if ans.upper() == "Y":
            self.saveInit()
        continueText()

    def saveInit(self):
        with open(f"./Save/{self.simName}.evsm", mode="w", newline="") as file:
            writer = c.writer(file)
            writer.writerow(["Parameter", "Value"])
            for key, value in self.parameters.items():
                writer.writerow([key, value])
        print(f"Parameters saved to ./Save/{self.simName}.evsm")

    def queryResults(self):
        clear()
        ans = input("Save results? (Y/N): ")
        while ans.upper() not in ["Y", "N"]:
            print("Please enter Y or N.")
            ans = input("Save initial parameters? (Y/N): ")

        if ans.upper() == "Y":
            self.saveResults()
        continueText()

    def saveResults(self):
        with open(f"./Results/results-{self.simName}-{d.today().strftime('%d-%m-%Y')}.txt", mode="w", newline="") as file:
            file.write("Simulation Results\n")
            for record in self.data:
                file.write(f"Generation {record['generation']}:\n")
                file.write(f"Population Size: {record['populationSize']}\n")
                file.write(f"Average Fitness: {record['averageFitness']:.2f}\n")
                file.write(f"Genotype Distribution: {record['genotypeDistribution']}\n")
                for trait, frequencies in record['phenotypeFrequency'].items():
                    file.write(f" {trait}: {frequencies}\n")
                file.write("\n")
        print(f"Results saved to ./Results/results-{self.simName}-{d.today().strftime('%d-%m-%Y')}.txt")
    
    def performGeneticAlgorithm(self):
        newPopulation = []
        populationSize = len(self.population)
    
        # Keep a percentage of the fittest entities
        eliteCount = int(populationSize * self.elitismRate)
        sortedPopulation = sorted(self.population, key=self.calculateFitness, reverse=True)
        newPopulation.extend(sortedPopulation[:eliteCount])
    
        # Fitness based crossover
        fitnessSum = sum(self.calculateFitness(genome) for genome in self.population)

        if fitnessSum == 0:
            weights = [1.0 / populationSize] * populationSize
        else:
            weights = [self.calculateFitness(genome) / fitnessSum for genome in self.population]
    
        for i in range(populationSize):
            parents = r.choices(self.population, weights=weights, k=2)
            parent1, parent2 = parents
            child = parent1.crossover(parent2)
            self.mutate(child)
            newPopulation.append(child)

        return newPopulation

    def mutate(self, genome):
        for genotype in genome.genotypes:
            newGeneticID = ""
            for allele in genotype.geneticID:
                if allele == "0" and r.random() < self.mutationRate:
                    newGeneticID += "1"
                elif allele == "1" and r.random() < self.mutationRate:
                    newGeneticID += "0"
                else:
                    newGeneticID += allele
            genotype.geneticID = genome.normaliseID(newGeneticID)

    def calculateFitness(self, genome):
        fitness = 0
        traitWeights = {"Speed": 1.0, "Metabolism": 1.2, "Bravado": 1.5}
        for i in range(len(genome.genotypes)):
            fitness += int(genome.genotypes[i].geneticID, 2) * traitWeights.get(self.phenotypes[i].trait, 1.0)
        return fitness

    def updatePopulation(self, newPopulation):
        self.population = newPopulation
        survivalThreshold = 2
        self.population = [genome for genome in self.population if self.calculateFitness(genome) > survivalThreshold]

    def collectData(self):
        populationSize = len(self.population)
        totalFitness = sum(self.calculateFitness(genome) for genome in self.population)
        averageFitness = totalFitness / populationSize
            
        genotypeDistribution = {}
        for genome in self.population:
            for genotype in genome.genotypes:
                if genotype.geneticID in genotypeDistribution:
                    genotypeDistribution[genotype.geneticID] += 1
                else:
                    genotypeDistribution[genotype.geneticID] = 1
            
        phenotypeFrequency = {phenotype.trait: {} for phenotype in self.phenotypes}
        for genome in self.population:
            for i in range(len(genome.genotypes)):
                genotype = genome.genotypes[i]
                phenotype = self.phenotypes[i]
                trait = genotype.getTrait()
                if trait not in phenotypeFrequency[phenotype.trait]:
                    phenotypeFrequency[phenotype.trait][trait] = 0

                phenotypeFrequency[phenotype.trait][trait] += 1

        self.data.append({
            "generation": self.generation,
            "populationSize": populationSize,
            "averageFitness": averageFitness,
            "genotypeDistribution": genotypeDistribution,
            "phenotypeFrequency": phenotypeFrequency
        })

    def displayStatistics(self):
        clear()
        latestData = self.data[-1]
        print(f"Generation {self.generation}:")
        print(f"Population Size: {latestData['populationSize']}")
        print(f"Average Fitness: {latestData['averageFitness']:.2f}")
        print(f"Genotype Distribution: {latestData['genotypeDistribution']}")
        print(f"Phenotype Frequency:")
        for trait, frequencies in latestData['phenotypeFrequency'].items():
            print(f" {trait}: {frequencies}")
        continueText()
        

class SimulationMenu:
    def __init__(self):
        self.running = True

    def displayMenu(self):
        while self.running:
            clear()
            print("EvoSim\n")
            print("1. Start Simulation")
            print("2. Load Simulation")
            print("3. Tutorial")
            print("4. Learn")
            print("5. Exit")
            choice = input(">>> ")
            self.handleChoice(choice)

    def handleChoice(self, choice: str):
        match choice:
            case "1":
                self.startSimulation()
            case "2":
                self.loadSimulation()
            case "3":
                self.showTutorials()
            case "4":
                self.showLearn()
            case "5":
                self.exitProgram()
            case _:
                print("Please chose a number between 1 and 5.")
                continueText()

    def startSimulation(self):
        self.setParameters()
        phenotypes = self.initialisePhenotypes()
        simulation = Simulation(self.mutationRate, phenotypes, self.elitismRate, self.initialPopulationSize, self.simName, self.genNum)
        self.displayInitialPopulation(simulation)
        simulation.runSimulation(self.genNum)


    def setParameters(self):
        clear()
        print("Set Simulation Parameters:")
        self.simName = input("Enter simulation name (e.g: example1): ")
        self.mutationRate = self.getFloatInput("Enter mutation rate (e.g: 0.01): ", 0, 1)
        self.elitismRate = self.getFloatInput("Enter elitism rate (e.g: 0.4): ", 0, 1)
        self.initialPopulationSize = self.getIntInput("Enter initial population size (e.g: 6): ", 2, 12)
        self.genNum = self.getIntInput("Enter number of generations (e.g: 4): ", 1, 9)
        continueText()
    
    def getFloatInput(self, prompt, minValue, maxValue):
        while True:
            try:
                value = float(input(prompt))
                if minValue is not None and (value < minValue or value > maxValue):
                    print(f"Please enter a float between {minValue} and {maxValue}.")
                else:
                    return value
            except ValueError:
                print(f"Please enter a float between {minValue} and {maxValue}.")
    
    def getIntInput(self, prompt, minValue, maxValue):
        while True:
            try:
                value = int(input(prompt))
                if minValue is not None and (value < minValue or value > maxValue):
                    print(f"Please enter a float between {minValue} and {maxValue}.")
                else:
                    return value
            except ValueError:
                print(f"Please enter a float between {minValue} and {maxValue}.")

    def initialisePhenotypes(self):
        speed = Phenotype("Speed", "s", ["Slow", "Fast"]) 
        metabolism = Phenotype("Metabolism", "m", ["Slow Metabolism", "Fast Metabolism"]) 
        bravado = Phenotype("Bravado", "b", ["Timid", "Courageous"])
        return [speed, metabolism, bravado]

    def displayInitialPopulation(self, simulation):
        clear()
        print("Initial Population:")
        for genome in simulation.population:
            genome.display()
            print()
        continueText()

    def loadParametersFromCSV(self, filename):
        parameters = {}
        with open(filename, mode="r") as file:
            reader = c.reader(file)
            next(reader)  # Skip the header
            for row in reader:
                parameters[row[0]] = row[1]
        return parameters

    def loadSimulation(self):
        clear()
        print("Load Simulation")
        
        # List all .evsm files in the Load folder
        files = [f for f in o.listdir("./Load") if f.endswith(".evsm")]
        
        if not files:
            print("No .evsm files found in the Load folder.")
            continueText()
            return

        print("Available .evsm files:")
        i = 1
        for file in files:
            print(f"{i}. {file}")
            i += 1
        
        choice = self.getIntInput("Enter the number of the file you want to load: ", 1, len(files))
        filename = files[choice - 1]
        filepath = f"./Load/{filename}"

        parameters = self.loadParametersFromCSV(filepath)
        
        self.simName = filename.split(".")[0]
        self.mutationRate = float(parameters["Mutation Rate"])
        self.elitismRate = float(parameters["Elitism Rate"])
        self.initialPopulationSize = int(parameters["Initial Population Size"])
        self.genNum = int(parameters["Number of Generations"])

        phenotypes = self.initialisePhenotypes()
        simulation = Simulation(self.mutationRate, phenotypes, self.elitismRate, self.initialPopulationSize, self.simName, self.genNum)
        self.displayInitialPopulation(simulation)
        simulation.runSimulation(self.genNum)

    def showTutorials(self):
        clear()
        print("Tutorial:")
        continueText()

    def showLearn(self):
        clear()
        print("Learn:")
        continueText()

    def exitProgram(self):
        clear()
        print("Exiting...")
        self.running = False

if __name__ == "__main__":
    createFolders()
    menu = SimulationMenu() 
    menu.displayMenu()
