import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from math import *
from scipy.integrate import odeint,solve_ivp
import matplotlib
from matplotlib.ticker import FuncFormatter, MultipleLocator
from matplotlib import rc


### Constantes physiques ###

g = 9.81  # Accélération due à la gravité

t0=0
dt=0.1
tf=10


### Paramètres initiaux sur le ballon

m = 0.450 # Masse du ballon en kg
a = 0.25  #longueur en m
b = 0.19  # largeur en m
Vb = 0.0048 #volume du ballon

alpha = 0  #angle d'attaque du coup de pied
gamma = 30 # angle de trajectoire de vole
khi = 0  #angle d'azimuth de la vitesse

w = 10 #norme de la vitesse de rotation en tour par seconde

  #theta wt à définir 

X0 = 20 #position initiale du ballon
Y0 = 50
Z0 = -0.5

U = 25  #composantes de la vitesse initiale
V = 1
W = 1
NV=sqrt(U**2+V**2+W**2) #norme de la vitesse

twt=np.arcsin(sqrt(V**2+W**2)/NV)

P = 10 #composante du vecteur vitesse angulaire (en s-1)
Q = 0
R = 0

PSI = 0
THET = 30
PHI = 0

Mij=np.array([[cos(THET)*cos(PSI), sin(THET)*cos(PSI)*sin(PHI)-cos(PHI)*sin(PSI), cos(PHI)*cos(PSI)*sin(THET)+sin(PHI)*sin(PSI)],[cos(THET)*sin(PSI), sin(THET)*sin(PSI)*sin(PHI)+cos(PHI)*cos(PSI), cos(PHI)*sin(PSI)*sin(THET)-sin(PHI)*cos(PSI)],[-sin(THET), sin(PHI)*cos(THET),cos(PHI)*cos(THET)]], dtype=float)



#### Definition des composantes (Xa, Ya, Za) de la force aéro ###

Cd = 2.44e-3-1.09e-3*twt+3.03e-4*twt**2+3.59e-7*twt**3-2.50e-8*twt**4  #coefficient de trainé
Cl = 6.25e-3*twt + 2.41e-4*twt**2 - 3.44e-6*twt**3#coefficient de portance 
Cm = 1.51e-2*twt - 1.69e-4*twt**2
Cy = 1

L = 0    #lift (masse volumique air à 1.2 kg/m3)
D = 0.5*1.2*Cd*NV*NV*Vb**(2/3) #drag
M = 0.5*1.2*Cm*NV*NV*Vb #pitching moment
Y = 0.5*1.2*Cy*NV*NV*Vb**(2/3) #side force

Xa = -D*U/NV
Ya = -D*V/NV + Y*W/(sqrt(V**2+W**2)) ## L=0 ici
Za = -D*W/NV - Y*V/sqrt(V**2+W**2)

## Composantes (La, Ma, Na) du moment aéro ##

La = 0  #L=N=0
Ma = M*W/sqrt(V**2+W**2)
Na = -M*V/sqrt(V**2+W**2)

### Definiton des fonctions pour les équa diff à résoudre ###



def FP():
    dPdt = La/0.0026   #moment inertie
    return dPdt

def FQ():
    dQdt = Ma/0.0033 + P*R*0.212121
    return dQdt

def FR():
    dRdt = Na/0.0033 - P*Q*0.212121
    return dRdt


def FPSI():
    dPSIdt = (Q*sin(PHI) + R*cos(PHI))/cos(THETA)
    return dPSIdt

def FTHET():
    dTHETdt = Q*cos(PHI) - R*sin(PHI)
    return dTHETdt

def FPHI():
    dPHIdt = P + tan(THET)*(Q*sin(PHI) + R*cos(PHI))
    

def FU():
    dUdt = Xa/m - g*sin(THET) - Q*W+ R*V
    return dUdt

def FV():
    dVdt = Ya/m + g*cos(THET)*sin(PHI) - R*U + P*W
    return dVdt

def FW():
    dWdt = Za/m + g*cos(THET)*cos(PHI) - P*V + Q*U


def Fpos(t,A,B,C):
    return np.dot(B,C)
    

### Initialisation des listes pour stocker les variables


vX=[X0]
vY=[Y0]
vZ=[Z0]

vU=[U]
vV=[V]
vW=[W]

vP=[P]
vQ=[Q]
vR=[R]

vPSI=[PSI]
vTHET=[THET]
vPHI=[PHI]

t_val=[t0]


### definition d'une intégration RK4 entre deux pas de temps

def runge_kutta(t, A, B, C, dt):
        k1 = dt * Fpos(t, A, B, C)
        k2 = dt * Fpos(t + 0.5*dt, A + 0.5*k1, B, C)
        k3 = dt * Fpos(t + 0.5*dt, A + 0.5*k2, B, C)
        k4 = dt * Fpos(t + dt, A + k3, B, C)
        A_new = A + (k1 + 2*k2 + 2*k3 + k4) / 6
        return A_new


## début de la boucle de résolution

while t_val[-1]<tf:
    
    t=t_val[-1]
    t_span=(t, t+dt)
    
    X=vX[-1]
    Y=vY[-1]
    Z=vZ[-1]
    
    U=vU[-1]
    V=vV[-1]
    W=vW[-1]
    
    P=vP[-1]
    Q=vQ[-1]
    R=vR[-1]
    
    PSI=vPSI[-1]
    THET=vTHET[-1]
    PHI=vPHI[-1]
    
    
    ####################### intégration ####################
    
    #On intègre notre position pour les vitesses et angles d'Euler initiaux déjà connu
    
    M_vit=np.array([[vU[-1]],[vV[-1]],[vW[-1]]])
    A=np.array([[vX[-1]],[vY[-1]],[vZ[-1]]])  #matrice des positions du ballon au temps t
    A_new=runge_kutta(t, A, Mij, M_vit, dt)   #matrice des positions du ballon au temps t+dt
    
    
    
    
    #update des valeurs de Xa, La, etc...    
    
    
    
    
    t_val.append(t+dt)
    
    
    
    
    





