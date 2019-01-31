import numpy as np
import matplotlib.pyplot as plt


class Pool:
    def __init__(self, summerTemp=23, winterTemp=36, tempThreshold=20, insulationThickness=0):
        self.summerTemp = summerTemp  # celsius
        self.winterTemp = winterTemp  # celsius
        self.summerDays = []
        self.summerTempThreshold = tempThreshold  # celsius
        self.insulationThickness = insulationThickness  # meters
        self.insulationPrice = 0
        self.surfaceArea = 200  # m^2
        self.sideArea = 90  # m^2
        self.windVelocity = 3  # m/s
        self.relativeHumidity = 0.5

        self.summerOff = True
        self.summerAirTemp = True

        self.nightCover = True
        self.nightCoverCost = 0
        self.coverReduction = 0.5 * 8/24
        self.standAlone = False

        self.hAirTop = 6.61  # W/(m^2 K)
        self.hAirWalls = 6.75  # W/(m^2 K)
        self.hmAir = 6.4133  # m/s
        self.vaporEnthalpy = 2418  # kJ/kg

        self.fig = None
        self.tempPlot = None
        self.heatPlot = None
        self.timeVector = np.linspace(0, 365, 365)
        self.heatLoss = np.zeros(365)  # W
        self.airTemp = None
        self.poolTemp = None
        self.evapMassRate = 0

        self.QEvap = None  # W ...
        self.QWaterInput = None
        self.QTop = None
        self.QWalls = None
        self.QRadiation = None

    def getTotalLoss(self):
        self.setAirTemp()
        self.setPoolTemp()

        self.getQEvap()
        self.getQWaterInput()
        self.getQSurfaces()
        self.getQRadiation()

        self.getCost()
        if self.standAlone:
            self.printStats()
            self.setPlot()

    def getLocalLoss(self):
        self.heatLoss = 0

        self.getQEvap()
        self.getQWaterInput()
        self.getQSurfaces()
        self.getQRadiation()

        return self.heatLoss

    def setAirTemp(self):
        self.airTemp = 6.4 + 29.5 * np.sin(1.5 * np.pi + 2 * np.pi * self.timeVector / 365)

    def setPoolTemp(self):
        self.summerDays = np.where(self.airTemp > self.summerTempThreshold)[0]
        if self.standAlone:
            print("Summer Pool for {} days a year.".format(len(self.summerDays)))

        self.poolTemp = np.ones(365) * self.winterTemp

        if self.summerAirTemp:
            self.poolTemp[self.summerDays] = self.airTemp[self.summerDays]
        else:
            self.poolTemp[self.summerDays] = self.summerTemp

    def getQEvap(self):
        surfaceSatVapor = self.getSatVaporAt(self.poolTemp)
        airSatVapor = self.getSatVaporAt(self.airTemp)
        airVapor = self.relativeHumidity * airSatVapor

        self.evapMassRate = self.hmAir * self.surfaceArea * (surfaceSatVapor - airVapor)  # kg/s

        if self.nightCover:
            self.evapMassRate *= 1 - self.coverReduction

        self.QEvap = self.evapMassRate * self.vaporEnthalpy
        self.heatLoss += self.QEvap

    def getSatVaporAt(self, temp):
        return 1.9747*10**(-5)*temp**2 + 1.3257*10**(-4)*temp + 3.9866*10**(-3)

    def getQWaterInput(self):
        coldWaterCp = 4.198  # kJ/(kg K) at 7 celsius
        self.QWaterInput = self.evapMassRate * coldWaterCp * (self.poolTemp - 7)
        self.heatLoss += self.QWaterInput

    def getQSurfaces(self):
        kGlass = 0.8  # W/mK
        kInsulation = 0.05  # W/mK

        resTop = 1 / (self.hAirTop * self.surfaceArea)
        if self.nightCover:
            resTop *= 1 + self.coverReduction

        resWalls = 0.015/kGlass + 1/self.hAirWalls
        if self.insulationThickness != 0:
            resWalls += self.insulationThickness/kInsulation

        resWalls /= self.sideArea
        resSum = resWalls + resTop

        resSurfaces = (1/resTop + 1/resWalls)**(-1)
        qTotal = (self.poolTemp - self.airTemp) / resSurfaces

        self.QTop = qTotal * (1 - resTop/resSum)
        self.QWalls = qTotal * (1 - resWalls/resSum)
        self.heatLoss += qTotal

    def getQRadiation(self):
        self.QRadiation = np.ones(365)
        self.QRadiation *= -8500

    def printStats(self):
        totalAnnualLoss = np.round(np.sum(self.heatLoss)*24/1000, 1)
        print("Total Heat Loss per year = {} kWh\n".format(totalAnnualLoss))

        for QLoss, label in zip([self.QEvap, self.QWaterInput, self.QTop, self.QWalls, self.QRadiation], ["Evap", "WaterInput", "Top surface", "Walls", "Radiation"]):
            annualLoss = np.round(np.sum(QLoss)*24/1000, 1)
            print("Q{} per year = {} kWh  ({}%)".format(label, annualLoss, np.round(annualLoss/totalAnnualLoss*100, 1)))

    def setPlot(self):
        self.fig, [self.tempPlot, self.heatPlot] = plt.subplots(2)
        self.tempPlot.datas, self.tempPlot.labels = [self.airTemp, self.poolTemp], ["Air", "Piscine"]
        self.heatPlot.datas, self.heatPlot.labels = [self.heatLoss, self.QEvap, self.QWaterInput, self.QTop, self.QWalls, self.QRadiation], ["Total", "Évaporation", "Remplissage", "Surface", "Murs", "Rayonnement"]

        for i, graph in enumerate([self.tempPlot, self.heatPlot]):
            for data, label in zip(graph.datas, graph.labels):
                if i == 0:
                    graph.set_title("Températures au cours de l'année")
                    graph.plot(self.timeVector, data, label=label)
                    graph.set_ylabel("T [$\degree$C]")
                else:
                    graph.set_title("Bilan des pertes thermiques annuelles")
                    graph.plot(self.timeVector, data/1000, label="{} ({}%)".format(label, np.round(100*np.sum(data)/np.sum(self.heatLoss), 1)))
                    graph.set_ylabel("$q$ [kW]")

            graph.set_xlabel("Jours")
            graph.legend(loc="best")
            graph.set_xlim(0, 365)

    def getCost(self):
        insulationVolume = self.sideArea*self.insulationThickness  # m^3
        self.insulationPrice = np.round(insulationVolume*100, 2)
        self.nightCoverCost = self.nightCover*50

        if self.standAlone:
            print("\nInsulation : {} m^3 => {} $".format(insulationVolume, self.insulationPrice))

    def showPlot(self):
        plt.tight_layout()
        plt.show()
