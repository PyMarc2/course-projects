from sympy import *
import numpy as np
import matplotlib.pyplot as plt

"""
Trouver la contrainte dans les cas suivant:
1. Configuration en vol de croisière (équilibre).
2. Configuration en équilibre avec un facteur de charge de 3g (poids x3).
3. Avion au sol soutenu par son train d’atterrissage.
"""


class AirPlane:
    def __init__(self):
        self.resolution = 2000

        self.wingLength = 150
        self.totalWeight = 0
        self.apparentWeight = 10
        self.nervures = np.array([[25, 50, 80.8, 111.5, 150], [7.13, 5.89, 6.86, 3.35, 3.21]])
        self.wheel = [25, 110] # cas 1 et 2
        # self.wheel = [25, -1815.93]
        self.bl = 1.648
        self.tl = 0.16

        self.chargeFactor = 1
        self.gravity = 1
        self.pAlum = 0.101
        self.pFuel = 0.0303
        self.beamAreaFunc = None
        self.V = []
        self.M = [0]

    def getWingShearAndMoment(self):
        # Fonction qui trace le graphique des efforts internes V(x) et M(x)
        self.getTotalWeight()
        self.getApparentWeight()

        A = np.linspace(150, 0, self.resolution)[1:]

        for i, a in enumerate(A):
            shear = self.getWingSectionWeight(a) - self.getAeroLoad(a) #AeroLoad = 0 lorsque l'avion est au sol (cas 3)
            self.V.append(shear)

        fig, (ax1, ax2) = plt.subplots(2)
        ax1.plot(A, self.V)
        ax1.set_ylabel("Effort tranchant [lb]")
        ax1.set_xlabel("Distance sur l'aile [po]")

        X, self.V = list(reversed(A)), list(reversed(self.V))
        Xstep = 150 / (self.resolution-1)

        for y in self.V:
            self.M.append((y*Xstep + self.M[-1]))

        for i in range(len(self.M)):
            self.M[i] += abs(min(self.M))

        ax2.plot(X, self.M[1:])
        ax2.set_ylabel("Moment de flexion [lb$\cdot$po]")

        plt.show()

    def getTotalWeight(self):
        # Fonction qui calcule le poids total de l'avion
        mainWeight = 5200 * self.gravity
        self.totalWeight = 2 * self.getWingSectionWeight() + mainWeight

    def getApparentWeight(self):
        # Fonction qui calcule le poids apparent (pour le cas 3)
        self.apparentWeight = self.totalWeight * self.chargeFactor

    def getWingSectionWeight(self, a=0, b=150):
        # Méthode des sections pour déterminer l'effort tranchant à la section voulue
        fuelWeight = self.getFuelVolume() * self.gravity * self.pFuel
        beamWeight = 2 * self.getBeamVolume() * self.gravity * self.pAlum
        coatWeight = self.getCoatingWeight(a, b)
        nervWeight = sum(self.nervures[1][np.where((self.nervures[0] >= a) & (self.nervures[0] <= b))[0]]) * self.gravity
        wheelWeight = self.wheel[1] * self.gravity if a <= self.wheel[0] else 0

        # print([fuelWeight, beamWeight, coatWeight, nervWeight, wheelWeight])

        weight = sum([fuelWeight, beamWeight, coatWeight, nervWeight, wheelWeight])

        return weight*self.chargeFactor

    def getFuelVolume(self, a=0, b=150):
        # Fonction qui calcule le volume de carburant pour la section voulue
        x = symbols('x')
        beamDistance = 14.4 - 4.8 * x * 10**-2
        beamHeight = 3.296 - 1.08967 * x * 10**-2
        beamAreaFunc = ((beamHeight-2*self.tl)*self.tl)+2*self.bl*self.tl

        areaFunc = beamDistance * beamHeight - 2*beamAreaFunc

        fuelVolume = integrate(areaFunc, (x, a, b))
        # meanX = integrate(x*areaFunc, (x, a, b)) / fuelVolume

        return fuelVolume  # , meanX

    def getBeamVolume(self, a=0, b=150):
        # Fonction qui caclule le volume des longerons à la section voulue
        x = symbols('x')
        beamHeight = 3.296 - 1.08967 * x * 10**-2
        self.beamAreaFunc = ((beamHeight-2*self.tl)*self.tl)+2*self.bl*self.tl

        beamVolume = integrate(self.beamAreaFunc, (x, a, b))
        # meanX = integrate(x*beamAreaFunc, (x, a, b)) / beamVolume

        return beamVolume  # , meanX

    def getAeroLoad(self, a=0, b=150):
        # Fonction qui calcule la force portance (intègre la fonction de lift donnée)
        x = symbols('x')
        liftFunc = 9 * self.apparentWeight / (16 * self.wingLength) * (1 - (x / self.wingLength)**8)
        aeroLoad = integrate(liftFunc, (x, a, b))
        # meanX = integrate(x*liftFunc, (x, a, b)) / aeroLoad

        return aeroLoad  # , meanX

    def getCoatingWeight(self, a=0, b=150):
        # Fonction qui calcule le poids du coating de l'aile
        x = symbols('x')
        weightFunc = 2.004 * 0.125 * self.pAlum * self.gravity * ((-20.6*x/150) + 41.2)

        coatWeight = integrate(weightFunc, (x, a, b))
        # meanX = integrate(x*weightFunc, (x, a, b)) / coatWeight

        return coatWeight  # , meanX

    def getInertia(self, x):
        # Fonction qui calcule l'inertie et le centroïde de l'aile
        bh = 3.296 - 1.08967 * x * 10 ** -2
        Y = bh/2
        X = (0.08*0.16*bh + 2*((1.648-0.16)/2+0.16)*1.488*0.16)/(0.16*bh + 2*(1.488*0.16))
        beamCentroid = (X, Y)
        print("Beam centroid:" + str(beamCentroid))

        

        # Iz = (((1.648*bh**3)/12) + (1.648*bh*((1.648/2) - X))**2) - ((1.488 * ((bh-0.32)**3)/12) + 1.488 * (bh-0.32)*((0.16 + 1.488/2)-X)**2)
        # Iz = 0.0018394*((-20.6*x/150) + 41.2)**3*0.125 + 2*(1.508834*10**(-3)*x + 0.452651 - 0.124*(-0.0109866*x+2.976)**3)
        Iz = (0.0018394 * 0.125 * (((-20.6/150) * x) + 41.2)**3) + 2*(((1/12) * self.bl * (bh**3)) - ((1/12) * (self.bl - self.tl) * (bh - 2*self.tl)**3))

        print("Inertia of beam @ x=%d :" % x + str(Iz))

        return Iz, X, Y

    def getNormalStress(self):
        # Fonction qui calcule et qui trace la contrainte normale
        X = np.linspace(0, 150, self.resolution)[1:]

        self.normalStress = []

        for i in range(len(X)):

            Iz, A, Y = self.getInertia(X[i])
            self.normalStress.append((-self.M[i] * Y / Iz)/1000)

        fig, ax1 = plt.subplots(1)
        ax1.plot(X, self.normalStress)
        ax1.set_ylabel("Contrainte normale [ksi]")
        ax1.set_xlabel("Distance sur l'aile [po]")
        plt.show()

    def getShearStress(self):
        # Fonction qui calcule et qui trace la contrainte de cisaillement
        X = np.linspace(0, 150, self.resolution)[1:]
        # Q = lambda x : (3.296 - 1.08967 * x * 10 ** -2)**2/8 * 0.16  + 0.16*1.488*((3.296 - 1.08967 * x * 10 ** -2)/2 - 0.08)
        Q = lambda x : 2*(2.414134*10**(-6)*x**2 - 2.75634*10**(-3)*x + 0.5905816)
        #  Q = ((self.bl * self.tl * (beamHeight/2 + self.tl/2)) + (self.tl * (beamHeight/2 - self.tl) * ((beamHeight/2 - self.tl/2)/2)) )* 2

        self.shearStress = []

        for i in range(len(X)):
            Iz, A, Y = self.getInertia(X[i])
            self.shearStress.append(((self.V[i] * Q(X[i])) / (Iz*0.16))/1000)

        fig, ax1 = plt.subplots(1)
        ax1.plot(X, self.shearStress)
        ax1.set_ylabel("Contrainte en cisaillement [ksi]")
        ax1.set_xlabel("Distance sur l'aile [po]")
        plt.show()

    def getFactor(self):
        # Fonction qui calcule le facteur de sécurité
        normalFactor = 60/max(list(map(abs, self.normalStress)))
        shearFactor = 25/max(list(map(abs, self.shearStress)))
        print(normalFactor)
        print(shearFactor)

plane = AirPlane()
plane.getWingShearAndMoment()
plane.getNormalStress()
plane.getShearStress()
plane.getFactor()

