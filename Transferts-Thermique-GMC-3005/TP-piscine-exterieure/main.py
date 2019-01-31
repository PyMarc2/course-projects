
def testPool():
    from heatLoss import Pool
    from geothermy import GeoSystem
    from heater import Heater

    bigPool = Pool(winterTemp=36, tempThreshold=20, insulationThickness=0.04)
    bigPool.standAlone = True
    bigPool.getTotalLoss()

    geoSystem = GeoSystem(bigPool, 1825, 1, 100, 2)
    geoSystem.tmiBlock = 275
    geoSystem.numberOfWell = 30
    geoSystem.wellDepth = 200

    heater = Heater(geoSystem, waterExArea=43, numberOfYears=5)
    heater.standAlone = True
    heater.getQExchanger()

    bigPool.showPlot()


def optimize():
    from optimizer import Optimizer
    Optimizer(numberOfYears=5).optimize()

if __name__ == '__main__':
    # testPool()
    optimize()
