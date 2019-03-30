
'''
Il s'agit de  trouver le facteur de sécurité des ailes d'un avion déjà construit.
On veut trouver les facteurs pour 3 cas particuliers:
- Accélération 3G
- Accélération 1G
- Équilibre au sol

Premièrement, on veut trouver le point apparent de l'avion afin de calculer le profil
de portance et faire un DCL complet de l'aile.

Deuxièmement,

Troisièmement,

'''

from sympy import *
import numpy as np
import matplotlib.pyplot as plt
init_printing()


class Airplane:
    def __init__(self):
        # Lengths are in 'in'
        self.wingLength = 150
        # self.wingEndWidth = 20.6
        # self.wingStartWidth = 41.2
        # self.wingCoatingWidth = 0.125

        # Beams parameters
        self.tl0 = 0.16
        self.bl0 = 1.648
        self.hl0 = 3.296

        self.t1f = 0.16
        self.blf = 1.648
        self.hlf = 1.648


        # Weights are in 'lbs'
        self.centerWeight = 5200
        self.wingWeight = 0
        self.wingFuelWeight = 0

        self.apparentWeight = 0

        self.aeroLoad = 0
        self.aeroLoadCentroid = 0

        self.G = 1
        self.chargeFactor = 1

        self.pCoating = 0.101  # lb/in3
        self.pFuel = 0.0303

        self.resolution = 1000

        self.shear = []

        self.wheelposition = [25, 110]
        self.nervures = np.array([[25, 50, 80.8, 111.5, 150], [7.13, 5.89, 6.86, 3.35, 3.21]])

    # ===== CAS1 - VOL À L'ÉQUILIBRE (1G) ===== #
    def getWingShearAndMoment(self):
        self.getTotalWeight()
        self.getTotalApparentWeight()

        A = np.linspace(150, 0, self.resolution)[1:]


        for i, a in enumerate(A):
            shear = self.getSectionWeight(a) - self.getAeroLoad(a)
            self.shear.append(shear)

        fig, (ax1, ax2) = plt.subplots(2)
        ax1.plot(A, self.shear)
        ax1.set_ylabel("Effort tranchant [lb]")
        ax1.set_xlabel("Distance sur l'aile [po]")

        X, self.shear = list(reversed(A)), list(reversed(self.shear))
        Xstep = 150 / (self.resolution-1)

        M = [0]
        for y in self.shear:
            M.append(y*Xstep + M[-1])

        ax2.plot(X, M[1:])
        ax2.set_ylabel("Moment")

        plt.show()

    def getSectionWeight(self, a=0):
        fuelWeight = self.getFuelWeight(a)
        beamWeight = self.getBeamWeight(a)
        coatWeight = self.getCoatingWeight(a)
        nervWeight = sum(
            self.nervures[1][np.where((self.nervures[0] >= a) & (self.nervures[0] <= self.wingLength))[0]]) * self.G
        wheelWeight = self.wheelposition[1] * self.G if a <= self.wheelposition[0] else 0

       # print([fuelWeight, beamWeight, coatWeight, nervWeight, wheelWeight])

        weight = sum([fuelWeight, beamWeight, coatWeight, nervWeight, wheelWeight])

        return weight

# === Poids apparent de l'avion === #
    def getAeroLoad(self, a=0):
        self.getTotalApparentWeight()
        x = symbols('x')
        w0 = self.apparentWeight
        L = self.wingLength

        self.weightFunctionAeroLoad = (9*w0/16*L)*(1-(x/L)**8)
        # print("Portance Répartition:", self.weightFunctionAeroLoad)

        self.totalAeroLoad = integrate(self.weightFunctionAeroLoad, (x, a, L))
        # print("Portance Totale:", self.totalAeroLoad)

        self.aeroLoadCentroid = integrate(self.weightFunctionAeroLoad*x, (x, a, L))/self.totalAeroLoad
        # print("Portance Centroid:", self.aeroLoadCentroid)

        return self.totalAeroLoad

    def getTotalApparentWeight(self):
        self.apparentWeight = self.getTotalWeight() * self.chargeFactor
        return self.apparentWeight

    def getTotalWeight(self):
        self.totalPlaneWeight = self.centerWeight + 2 * self.getWingWeight(a=0)
        print("Poids Total de l'avion:", self.totalPlaneWeight)
        return self.totalPlaneWeight

    def getWingWeight(self,a=0):
        self.getNervureWeight()
        self.wingTrainWeight = 110
        # print(self.wingTrainWeight)
        self.wingWeight = self.getCoatingWeight(a) + self.getBeamWeight(a) + self.getFuelWeight(a) + self.totalNervureWeight + self.wingTrainWeight
        # print("Poids d'une aile:", self.wingWeight)
        return self.wingWeight

    def getBeamWeight(self, a=0):
        x = symbols('x')
        beamHeight = 3.296 - 1.08967*x*10**-2
        updownVolume = self.bl0*self.tl0
        self.beamAreaFunc = (((beamHeight-2*self.tl0)*self.tl0)+2*updownVolume)
        self.totalBeamVolume = integrate(self.beamAreaFunc, (x, a, self.wingLength))

        self.weightFunctionBeam = 2 * self.beamAreaFunc * self.pCoating * self.G
        # print("Beam Répartition:", self.weightFunctionBeam)

        self.totalBeamWeight = 2 * self.totalBeamVolume * self.pCoating * self.G
        # print("Beam Poids Total:", self.totalBeamWeight)

        self.beamWeightCentroid = integrate(self.weightFunctionBeam*x, (x, a, self.wingLength))/self.totalBeamWeight
        # print("Beam Poids Centroid:", self.beamWeightCentroid)

        return self.totalBeamWeight

    def getFuelWeight(self, a=0):
        x = symbols('x')
        beamDistance = 14.4 - 4.8 * x * 10**-2
        beamHeight = 3.296 - 1.08967 * x * 10**-2
        self.rectangleVolumeFunc = beamDistance * beamHeight
        self.totalRectangleVolume = integrate(self.rectangleVolumeFunc, (x, a, self.wingLength))
        self.totalFuelVolume = self.totalRectangleVolume - 2*self.totalBeamVolume

        self.totalFuelWeight = self.totalFuelVolume * self.pFuel * self.G
        # print("Fuel Total Weight:", self.totalFuelWeight)

        self.weightFunctionFuel = (self.rectangleVolumeFunc - 2*self.beamAreaFunc) * self.pFuel * self.G
        # print("Fuel Répartition:", self.weightFunctionFuel)

        self.fuelWeightCentroid = (integrate(self.weightFunctionFuel*x, (x, a, self.wingLength)))/self.totalFuelWeight
        # print("Fuel Poids Centroid:", self.fuelWeightCentroid)

        return self.totalFuelWeight

    def getCoatingWeight(self, a=0):
        x = symbols('x')
        self.weightFunctionCoating = 2.004 * (0.125) * self.pCoating * ((-20.6*x/150) + 41.2)
        # print("Coating Répartition:", self.weightFunctionCoating)

        self.totalCoatingWeight = integrate(self.weightFunctionCoating, (x, a, self.wingLength))
        # print("Coating Poids Total:", self.totalCoatingWeight)

        self.coatingWeightCentroid = integrate(self.weightFunctionCoating*x, (x, a, self.wingLength))/self.totalCoatingWeight
        # print("Coating Weight Centroid:", self.coatingWeightCentroid)

        return self.totalCoatingWeight

    def getNervureWeight(self):
        self.totalNervureWeight = (7.13+5.89+6.86+3.35+3.21)
        # print("Nervure Total Weight:", self.totalNervureWeight)
        return self.totalNervureWeight



avion = Airplane()
avion.getWingShearAndMoment()