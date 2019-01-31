import numpy as np
import itertools


class Heater:
    def __init__(self, geoSystem, waterExArea=70, numberOfYears=5):
        self.standAlone = False
        self.skipGeothermy = False
        self.numberOfYears = numberOfYears
        self.geoSystem = geoSystem
        self.pool = geoSystem.pool
        self.waterExArea = waterExArea  # m^2
        self.vaporExArea = 0  # m^2  (init value...)
        self.efficacity = 0.75

        self.QTotal = np.zeros(365)  # W ...
        self.QWaterEx = np.zeros(365)
        self.QVaporEx = None
        self.QPac = None

        self.WPacTotal = 0
        self.geoWCost = 0
        self.pacCost = 0
        self.wellCost = 0
        self.waterCost = 0
        self.vaporCost = 0
        self.totalCost = 0

    def getQExchanger(self):
        self.getQWaterEx()
        self.geoSystem.heater = self

        if not self.skipGeothermy:
            self.checkGeothermy()
        self.correctYears()
        self.QTotal += self.QWaterEx
        self.getQVaporEx()
        self.getCost()
        if self.standAlone:
            self.setPlot()

    def getQWaterEx(self):
        heatCoef = 240  # W/(m^2 K)
        tempIn = 55  # Celsius
        NTU = 1.5
        Cmin = heatCoef * self.waterExArea / NTU  # kJ/(K s)
        self.QWaterEx = 0.75 * Cmin * (tempIn - self.pool.poolTemp)
        QOverflow = np.where(self.QWaterEx > self.pool.heatLoss)[0]
        self.QWaterEx[QOverflow] = self.pool.heatLoss[QOverflow]

        if self.pool.summerOff:
            self.QWaterEx[self.pool.summerDays] = 0

        self.QPac = self.QWaterEx[:365]/0.75

    def checkGeothermy(self):
        self.geoSystem.generateTemporalHeatMap()

        self.QPac = self.geoSystem.Qpac

        self.QWaterEx = 0.75 * self.QPac

        self.WPacTotal = np.sum(self.geoSystem.WPac) * 24

    def correctYears(self):
        if len(self.pool.timeVector) != self.numberOfYears * 365:
            self.pool.timeVector = np.linspace(0, self.numberOfYears*365, self.numberOfYears*365)
            self.pool.heatLoss = np.tile(self.pool.heatLoss, self.numberOfYears)
            self.pool.poolTemp = np.tile(self.pool.poolTemp, self.numberOfYears)

            newSummerDays = []
            for year in range(self.numberOfYears):
                localDays = np.add(np.array(self.pool.summerDays), 365*year).tolist()
                newSummerDays = list(itertools.chain(newSummerDays, localDays))
            self.pool.summerDays = newSummerDays

        self.QTotal = np.tile(self.QTotal, self.numberOfYears)

    def getQVaporEx(self):
        heatCoef = 340  # W/(m^2 K)
        tempIn = 115  # Celsius
        NTU = 1.39

        QDiff = np.max(self.pool.heatLoss - self.QWaterEx)
        Cmin = QDiff / (0.75 * (tempIn - np.max(self.pool.poolTemp)))  # kJ/(K s)
        self.vaporExArea = np.round(Cmin * NTU / heatCoef, 2)
        if self.vaporExArea < 0:
            self.vaporExArea = 0
            Cmin = 0

        self.QVaporEx = 0.75 * Cmin * (tempIn - self.pool.poolTemp)
        QOverflow = np.where(self.QVaporEx+self.QWaterEx > self.pool.heatLoss)[0]
        self.QVaporEx[QOverflow] = self.pool.heatLoss[QOverflow] - self.QWaterEx[QOverflow]

        if self.pool.summerOff:
            self.QVaporEx[self.pool.summerDays] = 0
        self.QTotal += self.QVaporEx

    def setPlot(self):
        self.pool.tempPlot.change_geometry(3, 1, 1)
        self.pool.heatPlot.change_geometry(3, 1, 2)
        heaterPlot = self.pool.fig.add_subplot(313)

        heaterPlot.datas, heaterPlot.labels = [self.QWaterEx, self.QVaporEx, self.QTotal], ["Water Exchanger", "Vapor Exchanger", "Total Exchanger"]
        for data, label in zip(heaterPlot.datas, heaterPlot.labels):
            heaterPlot.plot(self.pool.timeVector, data/1000, label=label)

        self.pool.tempPlot.legend(loc='best')
        heaterPlot.set_xlabel("Jours")
        heaterPlot.set_ylabel("$q$ [kW]")
        heaterPlot.legend(loc="best")
        heaterPlot.set_xlim(0, 364*self.numberOfYears)

    def getCost(self):
        waterExPrice = np.round(50*(self.waterExArea/90)**0.68, 2)
        vaporExStockPrice = np.round(50*(self.vaporExArea/90)**0.68, 2)
        vaporExRunPrice = np.round(4*(sum(self.QVaporEx)*24/1000000), 2)

        self.waterCost = np.round(waterExPrice, 2)
        self.vaporCost = np.round(vaporExStockPrice + vaporExRunPrice, 2)
        self.geoWCost = np.round(self.WPacTotal * 5 / 1000000, 2)
        self.pacCost = np.round(23*(np.max(self.QPac)/1000)**0.65, 2)
        self.wellCost = np.round(self.geoSystem.numberOfWell * self.geoSystem.wellDepth * 50 / 100)
        self.totalCost = np.round(self.waterCost + self.vaporCost + self.geoWCost + self.pacCost + self.wellCost + self.pool.insulationPrice + self.pool.nightCoverCost, 2)

        if self.standAlone:
            print("Water Exchanger : {} m^2 => {} $".format(np.round(self.waterExArea, 1), waterExPrice))
            print("Vapor Exchanger : {} m^2 => {} + {} = {} $".format(np.round(self.vaporExArea, 1), vaporExStockPrice, vaporExRunPrice, vaporExStockPrice+vaporExRunPrice))
            print("ECHANGERS : {} $".format(self.waterCost+self.vaporCost))
            print("GEO POWER : {} $".format(self.geoWCost))
            print("TOTAL COST = {} $".format(self.totalCost))

