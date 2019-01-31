import numpy as np
from heater import Heater
from geothermy import GeoSystem
from heatLoss import Pool


class Optimizer:
    def __init__(self, numberOfYears=5,
                 ins=(0.04, 0.041, 0.01),  # .|  0.04 m
                 exArea=(40, 60, 2),  # ......|  43 m^2
                 tmi=(273, 279, 0.5),  # .....|  275 K
                 wellNb=(30, 30.5, 3),  # ....|  30
                 depth=(200, 201, 10)):  # ...|  200 m

        self.numberOfYears = numberOfYears
        self.numberOfDays = numberOfYears * 365

        self.pools = []
        self.heaters = []
        self.geoSystems = []
        self.results = []

        self.insulationThicknesses = np.arange(ins[0], ins[1], ins[2])
        self.exchangerAreas = np.arange(exArea[0], exArea[1], exArea[2])
        self.tmiBlocks = np.arange(tmi[0], tmi[1], tmi[2])
        self.wellCounts = np.arange(wellNb[0], wellNb[1], wellNb[2])
        self.wellDepths = np.arange(depth[0], depth[1], depth[2])

        print("Number of System = ", np.prod([len(self.insulationThicknesses), len(self.exchangerAreas), len(self.tmiBlocks), len(self.wellCounts), len(self.wellDepths)]))

    def optimize(self):
        self.generatePools()
        self.generateGeothermy()
        self.generateHeaters()
        self.compute()
        self.listResults()

    def generatePools(self):
        for thickness in self.insulationThicknesses:
            pool = Pool(insulationThickness=thickness)
            pool.getTotalLoss()
            self.pools.append(pool)

    def generateGeothermy(self):
        for wellDepth in np.tile(self.wellDepths, len(self.exchangerAreas)):
            for wellCount in self.wellCounts:
                for tmiBlock in self.tmiBlocks:
                    for pool in self.pools:
                        geoSystem = GeoSystem(pool, self.numberOfDays, 5, 100, 2)
                        geoSystem.tmiBlock = tmiBlock
                        geoSystem.numberOfWell = wellCount
                        geoSystem.wellDepth = wellDepth
                        self.geoSystems.append(geoSystem)

    def generateHeaters(self):
        nbOfSystem = len(self.geoSystems) // len(self.exchangerAreas)

        for i, exArea in enumerate(self.exchangerAreas):
            for geoSystem in self.geoSystems[i*nbOfSystem:(i+1)*nbOfSystem]:
                heater = Heater(geoSystem, waterExArea=exArea, numberOfYears=self.numberOfYears)
                self.heaters.append(heater)

    def compute(self):
        print("Estimated time: {} minutes".format(len(self.heaters) * 7/15 // 60 + (len(self.heaters) * 7/15 % 60)/60))

        import time
        timeStart = time.time()

        for i, heater in enumerate(self.heaters):
            print("System {}/{}".format(i+1, len(self.heaters)))
            heater.getQExchanger()

        print("Timer = ", time.time() - timeStart)

    def listResults(self):
        params = ["INS", "WATER A", "VAPOR A", "TMI Bk", "INS $", "XWater $", "XVapor $", "Well NxL", "WELLS $", "PAC $", "Wpac $", "TOTAL $"]
        template = "|{" + "}|{".join("{}:>{}".format(i, len(param)+2) for i, param in enumerate(params)) + "}|"
        print(template.format(*params))

        for heater in self.heaters:
            self.results.append([heater.pool.insulationThickness, heater.waterExArea, heater.vaporExArea, np.round(heater.geoSystem.tmiBlock, 2),
                                 heater.pool.insulationPrice, heater.waterCost, heater.vaporCost,
                                 "{}x{}".format(heater.geoSystem.numberOfWell, heater.geoSystem.wellDepth),
                                 heater.wellCost, heater.pacCost, heater.geoWCost, heater.totalCost])

        for result in sorted(self.results, key=lambda x: x[-1]):
            print(template.format(*result))
