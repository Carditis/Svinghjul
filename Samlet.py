import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import stats
#import Model

""" Data array """
mursten1 = pd.read_csv("2m11_0950.csv",sep = ",")
mursten2 = pd.read_csv("2m11_0957.csv",sep = ",")
mursten3 = pd.read_csv("2m11_1000.csv",sep = ",")
mursten4 = pd.read_csv("2m11_1004.csv",sep = ",")


murstenArr = [mursten1, mursten2, mursten3, mursten4]


sving1 = pd.read_csv("2s11_1023.csv",sep = ",")
sving2 = pd.read_csv("2s11_1031.csv",sep = ",")
sving3 = pd.read_csv("2s11_1039.csv",sep = ",")
sving4 = pd.read_csv("2s11_1046.csv",sep = ",")

hjulArr = [sving1, sving2, sving3, sving4]

""" Dictionaries """
#Mursten
om = {}
om2 = {}
am = {}
am2 = {}
tm = {}
tm2 = {}
rm = {}
rm2 = {}
tausysm = {}
tausysm2 = {}
taufrikm = {}
taufrikm2 = {}

#Svinghjul
oh = {}
ah = {}
th = {}
rh = {}

""" Variables """
#Sten
a = 0.088 #m
b = 0.224 #m
m_sten = 2.224 #kg
I_sten = m_sten * 1/12 * (a**2 + b**2) #kg * m^2

#Stang
m_stang = 0.2244 #kg
r_stang = 0.005 #m
I_stang = 1/2 * m_stang * r_stang**2  #kg * m^2

#Trisse
m_trisse = 0.015 #kg
r_trisse = 0.01 #m
I_trisse = m_trisse * r_trisse**2 #kg * m^2
    
#Plexi
N_plexi = 3
m_plexi = 0.057 #kg
r_plexi = 0.144/2 #m
I_plexi = N_plexi * m_plexi * r_plexi**2 #kg * m^2

#Møtrik
m_møtrik = 0.01 #kg
r1_møtrik = 0.005 #m
r2_møtrik = 0.008425 #m
N_møtrik = 11
I_møtrik = 1/2 * m_møtrik * (r1_møtrik * r1_møtrik + r2_møtrik * r2_møtrik) * N_møtrik #kg * m^2

#Spændeskive
m_skive = 0.002 #kg
r1_skive = 0.005 #m
r2_skive = 0.009875 #m
N_skive = 9
I_skive = 1/2 * m_skive * (r1_skive * r1_skive + r2_skive * r2_skive) * N_skive #kg * m^2

#Svinghjul
m_sving = 4.980 #kg
I_sving = 0.5817162071 #kg * m^2

#Flaske
m_flaske = 2.078 #kg
g = 9.82 #N/kg
F_flaske = m_flaske * g #N
r_akse = 0.005 #m
tau_flaske = F_flaske * r_trisse

#Inertimoment
I_opstilling = I_stang + I_trisse + I_plexi + I_møtrik + I_skive #kg * m^2
I_system = I_opstilling + I_sten

#Normalkraft
m_opstilling = m_stang + m_trisse + m_plexi * N_plexi\
+ m_møtrik * N_møtrik + m_skive * N_skive + m_flaske #kg

m_stensystem = m_opstilling + m_sten #kg
m_svingsystem = m_opstilling + m_sving #kg

N_stensystem = m_stensystem * g #N
N_svingsystem = m_svingsystem * g #N


"""Functions"""

def murstensberegner (mArr):
    
    for i in range(len(mArr)):
        
        """Vinkelhastighed- og acceleration af mursten"""
        om["omegaM" + str(i+1)] = []
        om2["omegaM" + str(i+1)] = []
        am["alphaM" + str(i+1)] = []
        am2["alphaM" + str(i+1)] = []
        tm2["tidM" + str(i+1)] = []
        
        tm["tidM" + str(i+1)] = mArr[i]['time after start [s]'].tolist()
        rm["rotnumM" + str(i+1)] = mArr[i]['rotation number'].tolist()
        
        omegaberegner(tm["tidM" + str(i+1)], om["omegaM"+str(i+1)])
        alphaberegner(tm["tidM" + str(i+1)], om["omegaM"+str(i+1)], am["alphaM"+str(i+1)])
        
        """Moment af systemet"""
        tausysm["tausystemM" + str(i+1)] = []
        
        tausystemberegner(I_system, am["alphaM" + str(i+1)], tausysm["tausystemM" + str(i+1)])
        
        """Friktionsmoment"""
        taufrikm["taufrikM" + str(i+1)] = []
        taufrikm2["taufrikM" + str(i+1)] = []
        
        taufrikberegner(tausysm["tausystemM" + str(i+1)], tau_flaske, taufrikm["taufrikM" + str(i+1)])
        
        
        
        """Afskæring af data"""
        
        j = len(am["alphaM" + str(i+1)])-1
        while am["alphaM" + str(i+1)][j] < 0:
            j -= 1
            
               
        q = 0
        while (q+1) < j:
            
            om2["omegaM"+str(i+1)].append(om["omegaM"+str(i+1)][q])
            am2["alphaM" + str(i+1)].append(am["alphaM" + str(i+1)][q])
            tm2["tidM" + str(i+1)].append(tm["tidM1"][q])
            taufrikm2["taufrikM" + str(i+1)].append(taufrikm["taufrikM" + str(i+1)][q])
            q += 1
        
        """Plots"""
        
        #Hele turen
        hastighedsplot(tm["tidM" + str(i+1)],om["omegaM" + str(i+1)], 1)
        accelerationsplot(tm["tidM" + str(i+1)],am["alphaM" + str(i+1)], 2)
        friktionsmomentsplot(om["omegaM" + str(i+1)], taufrikm["taufrikM" + str(i+1)], 5)
        
        #Kun før
        hastighedsplot(tm2["tidM" + str(i+1)],om2["omegaM" + str(i+1)], 3)
        accelerationsplot(tm2["tidM" + str(i+1)],am2["alphaM" + str(i+1)], 4)
        friktionsmomentsplot(om2["omegaM" + str(i+1)], taufrikm2["taufrikM" + str(i+1)], 6)


def hjulberegner (hArr):

    
    for i in range(len(hArr)):
        """Vinkelhastighed- og acceleration af svinghjul"""
        oh["omegaS" + str(i+1)] = []
        ah["alphaS" + str(i+1)] = []
        
        th["tidS" + str(i+1)] = hArr[i]['time after start [s]'].tolist()
        rh["rotnumS" + str(i+1)] = hArr[i]['rotation number'].tolist()
        
        omegaberegner(th["tidS" + str(i+1)], oh["omegaS"+str(i+1)])
        alphaberegner(th["tidS" + str(i+1)], oh["omegaS"+str(i+1)], ah["alphaS"+str(i+1)])
        
        hastighedsplot(th["tidS" + str(i+1)],oh["omegaS" + str(i+1)], 7)
#        plt.plot(Model.tid_list,Model.omega_list, color="blue")
        
        accelerationsplot(th["tidS" + str(i+1)],ah["alphaS" + str(i+1)], 8)
    
        

def omegaberegner(tid, omega):
    for i in range(len(tid)):
        if i == 0:
            omega.append((2 * math.pi)/tid[i])
        else:
            omega.append((2 * math.pi)/(tid[i]-tid[i-1]))
            
            
def alphaberegner(tid, omega, alpha):
    for i in range(len(omega)):
        if i == 0 :
            alpha.append(omega[i]/tid[i])
        else:
            alpha.append((omega[i]-omega[i-1])/(tid[i]-tid[i-1])) 
            
def tausystemberegner(I_system, alpha, tausystem):
    for i in range(len(alpha)):
        tausystem.append(I_system*alpha[i])
        
def taufrikberegner(tausystem, tauflaske, taufrik):
   for i in range(len(tausystem)):
       taufrik.append(-tausystem[i]+tau_flaske)
        
            
def hastighedsplot (tid, omega, k):
    plt.figure(k).suptitle("Vinkelhastighed over for tid")
    plt.plot(tid,omega)
    plt.xlabel('Tid [s]')
    plt.ylabel('ω [rad/s]')
    plt.show
    
def accelerationsplot (tid, alpha, k):
    plt.figure(k).suptitle("Vinkelacceleration over for tid")
    plt.plot(tid,alpha)
    plt.xlabel('Tid [s]')
    plt.ylabel('α [rad/s²]')
    plt.show
    
def friktionsmomentsplot (vinkelhastighed, friktionsmoment, k):
    plt.figure(k).suptitle("Friktionsmoment over for vinkelhastighed")
    plt.plot(vinkelhastighed,friktionsmoment)
    plt.xlabel('ω [rad/s]')
    plt.ylabel('τ_frik [N*m]')
    plt.show
            
def tidappender (tid1, tid2):
    tid2.append(tid1)

murstensberegner(murstenArr)

hjulberegner(hjulArr)




            
       
            
            