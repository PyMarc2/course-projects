# Design Piscine Extérieure UL

### Considérations climatiques et  argent :



- Vents moyen par jour de 3m/s.
- Température extérieure: $T_{ext}(t) = 6.4+29.5\left[1.5\pi + \frac{2\pi t}{365} \right]$.
- Humidité relative moyenne: 50%.
- Quand la température extérieure dépasse 20°C, la piscine n'est pas chauffée.
- Sous 20°C de température extérieure, la piscine est chauffée à 36°C
- on cacul le ratio poképiece/CAD avec le prix moyen d'un MWh d'élextricité au québec. source: http://www.hydroquebec.com/residentiel/espace-clients/compte-et-facture/comprendre-facture/tarifs-residentiels-electricite/tarif-d.html
- Le ratio est d'environ 1 poképièce pour 10 CAD
- le prix d'une toile est d'environ 50 poképièce source:https://www.poolsuppliescanada.ca/covers/solar/solar-covers/solar-covers-18-x-36-ft-rectangle/

### Puits géothermiques

- Le liquide utilisé dans les puits géothermique est de l'eau salée à 50g de sel/kg d'eau (propriétés à 4°C)
  - $\rho = 1039.5kg/m^3$
  - $C_p = 3.8k\ J/Kg\cdot K$
  - $\mu = 1.55\times 10^{-3} kg/m\cdot s$
  -  $\text{Freezing point} \approx -5°C$
  - Corrélation 8.60
  - Vitesse nécessaire pour Reynolds turbulent = 1m/s
  - Vitesse posée: 3m/s
  - $h_i = 8550$

- La propagation du fluide caloporteur est modélisé afin de caluler sa température de sorti en fonction de sa temprature d'entrée et selon sa profoneur
- Nous avons posé le rayon du duct ( contenant les duex pipes et le coulis) à 0.1 m mais ce paramétre est un input dans la déclaration de la classe Ground.
  La fonction de conduction 1D a été utilisée afin de modéliser la variation temporelle de la température du sol. Cela permet de trouver la flux net de chaleur entre le liquide calorifique et le sol en tout temps et en tout point. 

Une approche matricielle a été utilisée avec l'équation de conduction discrétisée.

Les conditions limites (la première et dernière valeur de la matrice) sont respectivement le flux de chaleur selon Tmi et la température du premier noeud dans le sol, ainsi que la température au loin dans le sol (10°C). 

Ensuite, on veut modéliser la propagation du liquide dans l'échangeur géothermique. On pose donc un Tmi initial, et on compare avec la température du sol pour trouver un q et trouver la température du liquide caloporteur $\Delta l$ plus loin. Avec un pas assez petit, on devrait pas 

les tyaux sont en polythylene k =0.38  source -> http://ascpro0.ascweb.org/archives/cd/2010/paper/CPRT192002010.pdf

### Échangeurs

Optimisé par ludo, sleon les besoins en énergies


### PAC

$q_{PAC} = q_{sol} + W *COP$

$q_{PAC} = \dot{m}c_p(Tmo-Tmi)  + W \left( 1 + 5\frac{Tmi - Tmin}{12}\right)$

### Absorbtion de l'eau

On peut montrer le calcul de l'absorbtion de l'eau qui est quasi nulle au final

On veut donc ajouter une surface qui tirera profit de la radiation du soleil

### Fond noir

- On pose un 1.37kW/h à la surface de la terre provenant du soleil.

- Le déplacement du soleil dans le ciel induit un réduction de puissance d'un facteur 1/sqrt(2).

- L'ensoleillement d'une journée est de 12h pour une journée dégagée.

- Selon https://www.currentresults.com/Weather/Canada/Cities/sunshine-annual-average.php?fbclid=IwAR0D9XPR8s93ZbgqtoiGUvdh9TIV5AE6n5WbtxGRL_fuIAbkR_-jIk50iiI, L'ensoleillement moyen annuel est de 41%

https://www.wolframalpha.com/input/?i=0.1*((1.05*10%5E10)*(x-309))%5E(1%2F4)+*+(x-309)+%3D+750-(0.98*(5.67*10%5E-8)*x%5E4)

L'équation transcendante, reliant l'émission radiative du soleil à l'absorbtion du fond de la piscine, peinturé à une avec un certaine couleur, a été résolu afin de connaitre la température du fond de la piscine et de choisir les bons matériau pour ne pas dépasser 41°C, soit le seuil d'inconfort à la chaleur.

Durant l'hiver, la température de la piscine a été posé à 36°C.

Avec plusieurs itérations sur cette équation, on peut poser un T voulu à 40°C, la variable devient donc l'absorbitivité et l'émissivité de la surface. Le matériau choisit pour obtenir ces résultats est de la pierre ou du béton noir.

La résolution avec une seule itération, montre une température du fond à 40°C.

Le coefficient de convection naturel de l'eau à été calculé avec:

 https://www.wolframalpha.com/input/?i=0.1*((1.05*10%5E10)*(313-309))%5E(1%2F4)

Ce qui donne un q total de 180W/m^2 pour 200m^2 = 36kW.

Pour une journée ensoleillée, 432kWh d'énergie sera transmis par le fond vers l'eau.

En moyenne, on pourrait considérer une moyenne quotidienne qui vaudrait:

$E = 0,41*432 + 0,59*0,3*432 = 254kWh/jour$



### Toile durant la Nuit





### Limites d'Optimization

- Aire de l'échangeur à l'eau

  > - Range déterminé si $Q_{total}$ est totalement donnée par l'échangeur à l'eau.
  > - Range: [10, 85]$m^2$

- Aire de l'échangeur à vapeur

  > - Range déterminé si $Q_{total}$ est totalement donné par l'échangeur à vapeur.
  > - Est une conséquence directe de $T_{miBlock}$ et de la quantité de $Q$ complémentaire nécessaire
  > - 

- Grosseur de la PAC

  > - Dépend du Temperature Block

- Température de blocking de la PAC

  > [-2, 5]°C

- Profondeur des puits géothermique

  > [20 - 200]m par tranche de 5m

- Nombre de puits geothermique

  > [2 - 50]

- Isolation utilisée

  > [0 - 5 ]cm

- 



### Hypothèses

- La variation des propriétés de l'eau (densité, cp) en fonction de la température est négligée entre 10°C et -2°C