# TP - Transferts thermiques (GMC-3005)
# 2. Résolution numérique de l'équation de conduction en régime transitoire
# Hypothèses: direction radiale aucune résistance de contact entre le sol et l'échangeur

import numpy as np
import matplotlib.pyplot as plt
import math as m


class GeoSystem:
    def __init__(self, pool, numberOfDays, timeStep, numberOfMeters, distanceStep, summerOn=True):

        #  ======== INSTANCES ======

        self.pool = pool
        self.heater = None

        #  ======== USER INPUT VARIABLES ======

        self.numberOfDays = numberOfDays
        self.timeStep = timeStep
        self.numberOfTimeStep = int(self.numberOfDays / self.timeStep)
        self.summerOn = summerOn
        self.numberOfMeters = numberOfMeters
        self.distanceStep = distanceStep
        self.numberOfNodes = int(numberOfMeters/distanceStep)

        #  ======== PROFESSOR PARAMETERS ======

        self.groundConductivity = 2.5 # W/(m2⋅K)
        self.groundDiffusivity = 2 * 10 ** (-6) #
        self.groundrhocp = self.groundConductivity / self.groundDiffusivity

        self.T0 = 283 # température du sol et du fluide au temps 0

        self.tmiBlock = 275 # température seuil

        self.pipeInternalDiameter = 0.015
        self.pipeExternalDiameter = 0.020
        self.pipeMaterial = "Polyéthylène"
        self.pipeConductivity = 0.4

        self.ductRadius = 0.15
        self.ductShapeFactor = 21.9 * ((2 * self.ductRadius) / self.pipeExternalDiameter) ** (-0.38)

        self.coulisConductivity = 1.2

        #  ======== VARIABLE PARAMETERS ======

        self.numberOfWell = 20
        self.wellDepth = 200
        self.wellDistanceStep = 5

        self.calorificLiquid = {"Nature": ["Water"], "Density": [1040], "Cp": [3.8* 10 ** 3 ], "Speed": [3], "h": [8500]}
        # density : kg/m³
        # Cp : J/(kg·K)
        # h : W/(m2•K)
        self.chosenCalorificLiquid = 0
        self.totalLineicResistance = np.round((1 / 2) * ((1 / (m.pi * self.pipeInternalDiameter * self.calorificLiquid["h"][self.chosenCalorificLiquid])) + ((m.log(self.pipeExternalDiameter / self.pipeInternalDiameter)) / (2 * m.pi * self.pipeConductivity)))+ (1 / (self.ductShapeFactor * self.coulisConductivity)), 4) # (K·m)/W

        self.mdotTotal = m.pi * ((self.pipeInternalDiameter ** 2) / 4) * self.calorificLiquid["Speed"][self.chosenCalorificLiquid] * self.calorificLiquid["Density"][self.chosenCalorificLiquid] * self.numberOfWell # kg/s



        self.dailyTemperatureVariation = np.full(self.numberOfTimeStep, 0) # initialisation de liste


        #  ======== OTHER DECLARATIONS ======

        self.copList = [] # initialisation de liste
        self.dailyTmoList = [] # initialisation de liste [K]
        self.TmoList = [283] # K
        self.WPac = [] #  initialisation de liste [Watts]
        self.WoverFlowDays = [] # initialisation de liste [jour]

        self.TmiList = []
        self.timeStepSecond = self.timeStep * 3600 * 24 # seconds

        #  ========  INITIALIZATION OF THE MATRIX EQUATION SYSTEM ======

        self.liveTempMatrix = np.full((self.numberOfNodes, 1), (self.T0 * (self.groundrhocp / self.timeStepSecond)))
        self.coeffMatrix = np.zeros((self.numberOfNodes, self.numberOfNodes))
        self.motherMatrix = np.zeros((self.numberOfNodes, self.numberOfTimeStep)) # matrice des températures pour tous les noeuds. Chaque jours, la solution est stocké dans cette matrice


    # ====== MAIN FUNCTION ========

    def generateTemporalHeatMap(self):
        self.averageQpac() # formattage de liste
        self.Qpac = np.tile(self.Qpac, self.numberOfDays//365)# formattage de liste
        self.buildCoefficientsMatrix() # construction de la matrice des coefficients
        tempForNextTime = [283] # condition initial pour le premier jour, température de 283 K dans le sol

        for times in range(self.numberOfTimeStep -1): # boucle d'itération sur tous les temps

            if times == 0: # Résolution de l'équation et cacul des variables pour la première journée avec les conditions initiales
                self.TmiList.append(-((self.Qpac[times])/(self.mdotTotal * self.calorificLiquid['Cp'][self.chosenCalorificLiquid])) + 283) # calcul de Tmi

                if self.TmiList[-1] < self.tmiBlock: # bloquage de la température Tmi pour continuer à opérer la PAC à la température seuil
                    block = self.tmiBlock
                    self.TmiList[times] = self.tmiBlock
                    self.Qpac[times] = self.calculateQPACCorrected(block, self.TmoList[times])
                    self.copList.append(self.calculateCOP(block))
                    self.WPac.append(self.calculateWPAC(block, 4, self.Qpac[times], self.copList[times]))

                else:
                    self.copList.append(self.calculateCOP(self.TmiList[times]))
                    self.WPac.append(self.calculateWPAC(self.TmiList[times], 4, self.Qpac[times], self.copList[times]))

                self.dailyTmoList = (self.calculateTmoList(self.TmiList[times], self.wellDepth, self.wellDistanceStep, 283)) # calcul de Tmo pendant sa propagation dans les puits
                self.TmoList.append(self.dailyTmoList[-1])
                TmoMean = np.mean(self.dailyTmoList) # moyenne pondérée de Tmo
                self.dailyTemperatureVariation[times] = (self.TmiList[-1] - self.TmoList[-1])



                self.liveTempMatrix[0, 0] += ((TmoMean - 283)/self.totalLineicResistance) / (m.pi * self.ductRadius * self.distanceStep)  # Mise à jour des conditions limites pour la résolution du prochain jour
                self.liveTempMatrix[-1, -1] = self.T0 # mise à jour de la condition limite à l'infini
                tempForNextTime = self.solveTempForNextTime(self.liveTempMatrix) # solution du système d'équation
                self.storeTempData(tempForNextTime, times) # stockage des données den la journée

            else: # Résolution de l'équation et calcul des variables pour toutes les autres journées avec comme conditions initiales les paramètres de la veille

                self.TmiList.append(-((self.Qpac[times])/(self.mdotTotal*self.calorificLiquid['Cp'][self.chosenCalorificLiquid])) + self.dailyTmoList[-1])

                if self.TmiList[-1] < self.tmiBlock:
                        block = self.tmiBlock
                        self.TmiList[times] = self.tmiBlock
                        self.Qpac[times] = self.calculateQPACCorrected(block, self.TmoList[times])
                        self.copList.append(self.calculateCOP(block))
                        self.WPac.append(self.calculateWPAC(block, 4, self.Qpac[times], self.copList[times]))

                else:
                    self.copList.append(self.calculateCOP(self.TmiList[times]))
                    self.WPac.append(self.calculateWPAC(self.TmiList[times], 4, self.Qpac[times], self.copList[times]))


                self.dailyTmoList = (self.calculateTmoList(self.TmiList[times], self.wellDepth, self.wellDistanceStep, tempForNextTime[0][0]))
                self.TmoList.append(self.dailyTmoList[-1])
                TmoMean = np.mean(self.dailyTmoList)
                self.dailyTemperatureVariation[times] = (self.TmiList[-1] - self.TmoList[-1])



                self.liveTempMatrix = np.round(tempForNextTime * (self.groundrhocp / self.timeStepSecond), 4) #formattage avec le nombre de décimales voulu

                self.liveTempMatrix[0, 0] += ((TmoMean - tempForNextTime[0][0])/self.totalLineicResistance) / (m.pi * self.ductRadius * self.distanceStep)# Mise à jour des conditions limites pour la résolution du prochain jour
                self.liveTempMatrix[-1, -1] = self.T0 # mise à jour de la condition limite à l'infini
                tempForNextTime = self.solveTempForNextTime(self.liveTempMatrix)
                self.storeTempData(tempForNextTime, times)

        self.WPac = np.array(self.WPac) # reformattage des listes de données
        self.expandWpacQpacCop() # reformattage des listes de données


    # ===== CALCULATING FUNCTIONS  ====

    def averageQpac(self):
        self.Qpac = np.resize(self.heater.QPac, (int(self.heater.QPac.shape[0]/self.timeStep),self.timeStep))
        self.Qpac = np.average(self.Qpac,1)

    def expandWpacQpacCop(self):
        self.Qpac = np.repeat(self.Qpac,self.timeStep)
        self.WPac = np.repeat(self.WPac, self.timeStep)
        self.copList = np.repeat(self.copList, self.timeStep)

    def buildCoefficientsMatrix(self):
        firstRow = self.buildFirstRowCoeffMatrix() # création da la première ligne de la matrice des coefficients (pour le premier noeud)
        lastRow = self.buildLastRowCoeffMatrix() # création da la dernière ligne de la matrice des coefficients (pour le dernier noeud)
        for i in range(self.numberOfNodes): # boucle permettant de créer toute la matrice tridiagonale selon les équations du document complément
            for j in range(self.numberOfNodes):
                if i == j:  # Au noeud i
                    self.coeffMatrix[i, j] = np.round(self.groundrhocp / self.timeStepSecond + (2 * self.groundConductivity) / (self.distanceStep ** 2), 6)
                if i == j + 1:  # Au noeud i-1
                    self.coeffMatrix[i, j] = -np.round((self.groundConductivity * ((i + 1) * self.distanceStep - (self.distanceStep / 2))) / ((i + 1) * self.distanceStep * self.distanceStep ** 2), 6)
                if i == j - 1:  # Au noeud i+1
                    self.coeffMatrix[i, j] = -np.round((self.groundConductivity * ((i + 1) * self.distanceStep + (self.distanceStep / 2))) / ((i + 1) * self.distanceStep * self.distanceStep ** 2), 6)
        self.coeffMatrix[0, :] = firstRow
        self.coeffMatrix[-1, :] = lastRow
        return np.round(self.coeffMatrix, 6)

    def buildFirstRowCoeffMatrix(self):
        firstCoeff = [(self.groundrhocp / self.timeStepSecond) + (2 * self.groundConductivity * (self.ductRadius + (self.distanceStep / 2))) /(self.ductRadius * self.distanceStep ** 2),- (2 * self.groundConductivity * (self.ductRadius + (self.distanceStep / 2))) / (self.ductRadius * self.distanceStep ** 2)] # selon l'équation du document de complément pour la première ligne
        listOfZeros = [0 for i in range(self.numberOfNodes - 2)] # formattage
        return np.round(np.concatenate((firstCoeff, listOfZeros)), 6)

    def buildLastRowCoeffMatrix(self):
        lastCoeff = [1] # condiion limites (pour obtenir 283 au dernier noeud)
        listOfZeros = [0 for i in range(self.numberOfNodes - 1)] # formattage
        return np.concatenate((listOfZeros, lastCoeff))

    def solveTempForNextTime(self, liveTempMatrix):
        tempNextTimeMatrix = np.linalg.solve(self.coeffMatrix, liveTempMatrix)
        return np.round(tempNextTimeMatrix, 5)

    def storeTempData(self, tempMatrix, absoluteTime):
        self.motherMatrix[:, absoluteTime] = tempMatrix[:, 0]

    def calculateTmoList(self, Tmi, wellDepth, depthStep, Tsol):
        dailyTmoList = [0 for i in range(int(2*wellDepth / depthStep))]
        for i in range(int(2*wellDepth / depthStep)):
            dailyTmoList[i] = np.round(((np.round(Tsol, 6)-np.round(Tmi, 6))/(self.totalLineicResistance * (self.mdotTotal/self.numberOfWell) * self.calorificLiquid["Cp"][self.chosenCalorificLiquid]))*depthStep + Tmi, 6)

            if Tsol-Tmi <= 10**(-4): # mise à zéro pour éviter les erreurs numériques
                Tmi = Tsol
            if dailyTmoList[i] >= Tsol: # bloquage de la température pour empêché de dépassement numériques
                dailyTmoList[i] = Tsol
                Tmi = dailyTmoList[i]
            else:
                Tmi = dailyTmoList[i]
        return dailyTmoList

    def calculateCOP(self, Tmi):
        cop = 1 + 5*((Tmi-271)/12)
        return cop

    def calculateWPAC(self, Tmi, Tmo, Qpac, cop):
        W = Qpac/cop
        return W

    def calculateQPACCorrected(self, Tmi, Tmo):
        Qpac = (self.mdotTotal * self.calorificLiquid['Cp'][self.chosenCalorificLiquid] * (Tmo - Tmi))
        return Qpac



     #===== GRAPHIC DISPLAY FUNCTIONS  ========

    def graphOneNodeAtAllTime(self, x):
        plt.plot(range(self.numberOfTimeStep), self.motherMatrix[x, :])
        plt.title("Temperature of node #{} for all {} days.".format(x, self.numberOfDays))
        plt.xlabel("[Days]")
        plt.ylabel("Temperature [°K]")
        plt.show()

    def graphTmoTmiAtAllTime(self):
        plt.plot(range(len(self.TmoList)), self.TmoList, 'k', label="Tmo")
        plt.plot(range(len(self.TmiList)), self.TmiList, 'r', label="Tmi")
        plt.title("Temperature of Calorific Liquid for all days")
        plt.legend()
        plt.xlabel("[Days]")
        plt.ylabel("Temperature [°K]")
        plt.show()

    def graphQSystems(self):
        plt.plot(range(len(self.Qpac)), self.heater.QVaporEx/1000, 'k', label="QVapor")
        plt.plot(range(len(self.Qpac)), self.heater.QWaterEx/1000, 'r', label="QWater")
        # plt.plot(range(len(self.Qpac)), self.heater.QPac/1000 - self.Qpac/1000, 'b', label="QVapor")
        plt.legend(loc=1)
        plt.title('Q of systems')
        plt.xlabel("[Days]")
        plt.ylabel("Heat [kW]")
        plt.show()

    def graphSuperpositionProfiles(self):
        for i in range(self.numberOfTimeStep-10):
            plt.plot(range(self.numberOfNodes), self.motherMatrix[:, i])
        plt.title("Profile Evolution")
        plt.xlabel("[Meters]")
        plt.ylabel("Temperature [°K]")
        plt.show()

    def graphWPACAllTime(self):
        fig, ax1 = plt.subplots()
        ax1.set_xlabel("Days")
        ax1.plot(range(len(self.WPac)), self.WPac/1000, 'b', label="WPac")
        ax2 = ax1.twinx()
        ax2.plot(range(len(self.copList)), self.copList, 'r', label="COP")
        ax1.legend(loc=1)
        ax2.legend(loc=2)
        plt.title("PAC Informations")
        ax2.set_ylabel("COP")
        ax1.set_ylabel("Wpac [kW]")
        plt.show()

    def graphSuperpositionNodes(self):
        for i in range(self.numberOfNodes):
            plt.plot(range(self.numberOfTimeStep - 10), self.motherMatrix[i, 0:self.numberOfTimeStep - 10])
        plt.title("Nodes Temperature Variation")
        plt.xlabel("[Days]")
        plt.ylabel("Temperature [°K]")
        plt.show()

    def printSimulationVariableParameters(self):
        print("\n ==================================================")
        print("(|              VARiABLE PARAMETERS               |)")
        print(" ==================================================")
        print("Simulation Time: {} Days".format(self.numberOfDays))
        print("Step Time: {} Days/Step".format(self.timeStep))
        print("Simulation Distance: {} Meters".format(self.numberOfMeters))
        print("Step Distance: {} Meters/Step".format(self.distanceStep))
        print("Well Depth: {} Meters".format(self.wellDepth))
        print("Number Of Wells: {} ".format(self.numberOfWell))
        print("Well Resolution for Tmo Calculation: {} Meters".format(self.wellDistanceStep))
        print("Calorific Liquid Nature: {}".format(self.calorificLiquid["Nature"][self.chosenCalorificLiquid]))
        print("Calorific Liquid Speed: {} m/s".format(self.calorificLiquid["Speed"][self.chosenCalorificLiquid]))
        print("Calorific Liquid Density: {} kg/m^3".format(self.calorificLiquid["Density"][self.chosenCalorificLiquid]))
        print("Calorific Liquid Specific Heat: {} J/kg*K".format(self.calorificLiquid["Cp"][self.chosenCalorificLiquid]))
        print("Calorific Liquid Convection Coefficient: {} W/m^2*K".format(self.calorificLiquid["h"][self.chosenCalorificLiquid]))
        print("Calorific Liquid Mass Transfer per Well: {} kg/s*Well".format(self.mdotTotal/self.numberOfWell))
        print("Total Calorific Liquid Mass Transfer: {} kg/s".format(self.mdotTotal))
        print("Summer Days: {}".format(self.summerOn))
        print("Total Duct Lineic Resistance: {} K*m/W".format(self.totalLineicResistance))

    def printSimulationInfos(self):
        print("\n ==================================================")
        print("(|            SiMULATiON iNFORMATiONS             |)")
        print(" ==================================================")
        print("Calorific Liquid Lowest Hot : {} °K".format(min(self.TmoList)))
        print("Calorific Liquid Lowest Cold : {} °K".format(min(self.TmiList)))
        print("Ground 1st Node Lowest Temperature: {} °K".format(min(self.motherMatrix[0, :])))


if __name__ == '__main__':
    from heatLoss import Pool
    from heater import Heater

    def testPool():
        bigPool = Pool(winterTemp=36, tempThreshold=20, insulationThickness=0.04)
        bigPool.getTotalLoss()

        geoSystem = GeoSystem(bigPool, 1825, 1, 100, 2)
        geoSystem.tmiBlock = 275
        geoSystem.numberOfWell = 30
        geoSystem.wellDepth = 200

        heater = Heater(geoSystem, waterExArea=43, numberOfYears=5)
        heater.getQExchanger()

        geoSystem.graphTmoTmiAtAllTime()
        geoSystem.graphQSystems()
        geoSystem.graphSuperpositionNodes()
        geoSystem.graphWPACAllTime()
        geoSystem.graphSuperpositionProfiles()

    testPool()


