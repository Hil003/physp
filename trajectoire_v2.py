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

alpha = 60  #angle d'attaque du coup de pied
gamma = 30 # angle de trajectoire de vol
khi = 0  #angle d'azimuth de la vitesse

w = 20 #norme de la vitesse de rotation en tour par seconde


X0 = 20 #position initiale du ballon
Y0 = 50
Z0 = 0

U = 32 #composantes de la vitesse initiale
V = 0.1
W = 0.1


P = 10 #composante du vecteur vitesse angulaire (en s-1)
Q = 0
R = 0

PSI = 0 #angles d'euler initiaux
THET = 20
PHI = 0



### Definiton des fonctions pour les équa diff à résoudre ###


def FP(La):
    dPdt = La/0.0026   #moment inertie
    return dPdt

def FQ(Ma,P,R):
    dQdt = Ma/0.0033 + P*R*0.212121
    return dQdt

def FR(Na,P,Q):
    dRdt = Na/0.0033 - P*Q*0.212121
    return dRdt


def FPSI(PHI, THET, Q, R):
    dPSIdt = (Q*sin(PHI) + R*cos(PHI))/cos(THET)
    return dPSIdt

def FTHET(PHI, Q, R):
    dTHETdt = Q*cos(PHI) - R*sin(PHI)
    return dTHETdt

def FPHI(THET, PHI, P, Q, R):
    dPHIdt = P + tan(THET)*(Q*sin(PHI) + R*cos(PHI))
    return dPHIdt
    
def FU(Xa,THET,Q,R,V,W):
    dUdt = Xa/m - g*sin(THET) - Q*W+ R*V
    return dUdt


def FV(Ya,THET,PHI,P,R,U,W):
    dVdt = Ya/m + g*cos(THET)*sin(PHI) - R*U + P*W
    return dVdt

def FW(Za,THET,PHI,P,Q,U,V):
    dWdt = Za/m + g*cos(THET)*cos(PHI) - P*V + Q*U
    return dWdt
    

def FX(Mij,U,V,W):
    C=np.array([[U],[V],[W]])
    return np.dot(Mij[0],C)

def FY(Mij,U,V,W):
    C=np.array([[U],[V],[W]])
    return np.dot(Mij[1],C)

def FZ(Mij,U,V,W):
    C=np.array([[U],[V],[W]])
    return np.dot(Mij[2],C)
    


### Résolution du système d'équation différentielle #####


def système(t, variables):
    
    X, Y, Z, U, V, W, P, Q, R, PSI, PHI, THET = variables
    
    #calcul les valeurs intermediaires 
    
    NV=sqrt(U**2+V**2+W**2) #norme de la vitesse

    twt=np.arcsin(sqrt(V**2+W**2)/NV)
    Cd = 2.44e-3-1.09e-3*twt+3.03e-4*twt**2+3.59e-7*twt**3-2.50e-8*twt**4  #coefficient de trainé
    Cl = 6.25e-3*twt + 2.41e-4*twt**2 - 3.44e-6*twt**3#coefficient de portance 
    Cm = 1.51e-2*twt - 1.69e-4*twt**2
    Cy = twt*(-1.50e-3 + 6.49e-4*w -8.35e-5*w**2) + (-4.11e-5 - 3.82e-5*w + 2.64e-6*w**2)*twt**2 + (4.94e-7 + 2.74e-7*w -1.77e-8*w**2)*twt**3
    
    
    L = 0    #lift (masse volumique air à 1.2 kg/m3)
    D = 0.5*1.2*Cd*NV*NV*Vb**(2/3) #drag
    M = 0.5*1.2*Cm*NV*NV*Vb #pitching moment
    Y = 0.5*1.2*Cy*NV*NV*Vb**(2/3) #side force
    
    
    Xa = -D*U/NV
    Ya = -D*V/NV + Y*W/(sqrt(V**2+W**2)) ## L=0 ici
    Za = -D*W/NV - Y*V/sqrt(V**2+W**2)
    
    La = 0  #L=N=0
    Ma = M*W/sqrt(V**2+W**2)
    Na = -M*V/sqrt(V**2+W**2)
    
    Mij=np.array([[cos(THET)*cos(PSI), sin(THET)*cos(PSI)*sin(PHI)-cos(PHI)*sin(PSI), cos(PHI)*cos(PSI)*sin(THET)+sin(PHI)*sin(PSI)],[cos(THET)*sin(PSI), sin(THET)*sin(PSI)*sin(PHI)+cos(PHI)*cos(PSI), cos(PHI)*sin(PSI)*sin(THET)-sin(PHI)*cos(PSI)],[-sin(THET), sin(PHI)*cos(THET),cos(PHI)*cos(THET)]], dtype=float)

    
    ## definition du système d'equa diff
    
    dX_dt=FX(Mij,U,V,W)
    dY_dt=FY(Mij,U,V,W)
    dZ_dt=FZ(Mij,U,V,W)
    
    dU_dt=FU(Xa,THET,Q,R,V,W)
    dV_dt=FV(Ya,THET,PHI,P,R,U,W)
    dW_dt=FW(Za,THET,PHI,P,Q,U,V)
    
    dP_dt=FP(La)
    dQ_dt=FQ(Ma,P,R)
    dR_dt=FR(Na,P,Q)
    
    dPSI_dt=FPSI(PHI, THET, Q, R)
    dPHI_dt=FPHI(THET, PHI, P, Q, R)
    dTHET_dt=FTHET(PHI, Q, R)
    

    return [dX_dt,dY_dt,dZ_dt,dU_dt,dV_dt,dW_dt,dP_dt,dQ_dt,dR_dt,dPSI_dt,dPHI_dt,dTHET_dt]



initial_variables=[X0,Y0,Z0,U,V,W,P,Q,R,PSI,PHI,THET]
num_points=100
t_eval = np.linspace(t0, tf, num_points)


solution=solve_ivp(système, (t0,tf), initial_variables, method='DOP853',t_eval=t_eval)


t = solution.t

X = solution.y[0]
Y = solution.y[1]
Z = solution.y[2]
U = solution.y[3]
V = solution.y[4]
W = solution.y[5]
P = solution.y[6]
Q = solution.y[7]
R = solution.y[8]
PSI = solution.y[9]
PHI = solution.y[10]
THET = solution.y[11]



plt.plot(X,-Z)
plt.ylim(bottom=0)

