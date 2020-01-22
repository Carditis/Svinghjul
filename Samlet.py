import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score




""" Data array """
mursten1 = pd.read_csv("2m11_0950.csv",sep = ",")
mursten2 = pd.read_csv("2m11_0957.csv",sep = ",")
mursten3 = pd.read_csv("2m11_1000.csv",sep = ",")
mursten4 = pd.read_csv("2m11_1004.csv",sep = ",")

mursten1 = pd.read_csv("1m11_1423.csv",sep = ",")

murstenArr = [mursten2, mursten3, mursten4]
# murstenArr = [mursten1]


# sving1 = pd.read_csv("2s11_1023.csv",sep = ",")
sving2 = pd.read_csv("2s11_1031.csv",sep = ",")
sving3 = pd.read_csv("2s11_1039.csv",sep = ",")
sving4 = pd.read_csv("2s11_1046.csv",sep = ",")

# sving1 = pd.read_csv("1s11_1402.csv",sep = ",")

hjulArr = [sving2, sving3, sving4]
# hjulArr = [sving1]

""" Dictionaries """
#Mursten
om = {}
om2 = {}
om3 = {}
am = {}
am2 = {}
am3 = {}
tm = {}
tm2 = {}
tm3 = {}
rm = {}
rm2 = {}
tausysm = {}
tausysm2 = {}
taufrikm = {}
taufrikm2 = {}
taufrikm3 = {}

#Svinghjul
oh = {}
oh2 = {}
oh3 = {}
ah = {}
ah2 = {}
ah3 = {}
th = {}
th2 = {}
th3 = {}
rh = {}
rh2 = {}
tausys = {}
svinginerti = {}



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
r_sving = 0.4 #m

#Flaske
m_flaske = 2.078 #kg
g = 9.82 #N/kg
F_flaske = m_flaske * g #N
r_akse = 0.005 #m
tau_flaske = F_flaske * r_trisse

#Inertimoment
I_opstilling = I_stang + I_trisse + I_plexi + I_møtrik + I_skive #kg * m^2
I_systemM = I_opstilling + I_sten


#Normalkraft
m_opstilling = m_stang + m_trisse + m_plexi * N_plexi\
+ m_møtrik * N_møtrik + m_skive * N_skive + m_flaske #kg

m_stensystem = m_opstilling + m_sten #kg
m_svingsystem = m_opstilling + m_sving #kg

N_stensystem = m_stensystem * g #N
N_svingsystem = m_svingsystem * g #N



"""Functions"""

def func(x1,a,b):
    return a*x1+b

def MurstensModel():
   
    modelM_omega = []
    modelM_omega.append(0)
    
    modelM_alpha = []
    modelM_alpha.append(0)
    
    modelM_tid = []
    modelM_tid.append(0)
    
    
    
    
    i = 1
    S = 0
    for j in range(int(tm["tidM1"][tilstandsskifteM]/(1/10))):
        modelM_tid.append(i/10)
        modelM_omega.append(modelM_omega[i-1] + modelM_alpha[i-1] * (1/10))
        modelM_alpha.append((tau_flaske - (slope_før * modelM_omega[i] + intercept_før))/I_systemM)
        deltaS = modelM_omega[i] * (1/10) * r_trisse
        S += deltaS
        i += 1
    
    while modelM_omega[i-1] > 0:
        modelM_tid.append(i/10)
        modelM_omega.append(modelM_omega[i-1] + modelM_alpha[i-1] * (1/10))
        modelM_alpha.append((-(slope_før * modelM_omega[i] + intercept_før))/I_systemM)
        i += 1
        
    hastighedsplot(modelM_tid, modelM_omega, 4)
    accelerationsplot(modelM_tid, modelM_alpha, 5)
    
    
# print(modelM_tid)
# print(modelM_omega)
    
def SvinghjulsModel():
   
    modelS_omega = []
    modelS_omega.append(0)
    
    modelS_alpha = []
    modelS_alpha.append(0)
    
    modelS_tid = []
    modelS_tid.append(0)
    
    I_hjulTeori = 0.5817162071 #N*m

    
    I_systemS = I_opstilling + I_hjulTeori
    
    i = 1
    S = 0
    for j in range(int(th["tidS1"][tilstandsskifteS]/(1/10))):
        modelS_tid.append(i/10)
        modelS_omega.append(modelS_omega[i-1] + modelS_alpha[i-1] * (1/10))
        modelS_alpha.append((tau_flaske - (slope_før * modelS_omega[i] + intercept_før))/I_systemS)
        deltaS = modelS_omega[i] * (1/10) * r_trisse
        S += deltaS
        i += 1
    
    while modelS_omega[i-1] > 0:
        modelS_tid.append(i/10)
        modelS_omega.append(modelS_omega[i-1] + modelS_alpha[i-1] * (1/10))
        modelS_alpha.append((-(slope_før * modelS_omega[i] + intercept_før))/I_systemS)
        i += 1
        
    hastighedsplot(modelS_tid, modelS_omega, 9)
    accelerationsplot(modelS_tid, modelS_alpha, 10)
    
    
# print(modelM_tid)
# print(modelM_omega)

def murstensberegner (mArr):
    
    merge_omega_list_før = []
    merge_taufrik_list_før = []
    merge_omega_list_efter = []
    merge_taufrik_list_efter = []

    for i in range(len(mArr)):

        """Vinkelhastighed- og acceleration af mursten"""
        om["omegaM" + str(i+1)] = []
        om2["omegaM" + str(i+1)] = []
        om3["omegaM" + str(i+1)] = []
        am["alphaM" + str(i+1)] = []
        am2["alphaM" + str(i+1)] = []
        am3["alphaM" + str(i+1)] = []

        tm["tidM" + str(i+1)] = mArr[i]['time after start [s]'].tolist()
        tm2["tidM" + str(i+1)] = []
        tm3["tidM" + str(i+1)] = []
        rm["rotnumM" + str(i+1)] = mArr[i]['rotation number'].tolist()


        omegaberegner(tm["tidM" + str(i+1)], om["omegaM"+str(i+1)])
        alphaberegner(tm["tidM" + str(i+1)], om["omegaM"+str(i+1)], am["alphaM"+str(i+1)])

        """Moment af systemet"""
        tausysm["tausystemM" + str(i+1)] = []

        tausystemberegner(I_systemM, am["alphaM" + str(i+1)], tausysm["tausystemM" + str(i+1)])

        """Friktionsmoment"""
        taufrikm["taufrikM" + str(i+1)] = []
        taufrikm2["taufrikM" + str(i+1)] = []
        taufrikm3["taufrikM" + str(i+1)] = []

        taufrikberegner(tausysm["tausystemM" + str(i+1)], tau_flaske, taufrikm["taufrikM" + str(i+1)])



        """Afskæring af data"""

        j = len(am["alphaM" + str(i+1)])-2
        while am["alphaM" + str(i+1)][j+1] < 0:
            j -= 1
        global tilstandsskifteM
        tilstandsskifteM = j


        q = 0
        while (q+1) < j:

            om2["omegaM"+str(i+1)].append(om["omegaM"+str(i+1)][q])
            am2["alphaM" + str(i+1)].append(am["alphaM" + str(i+1)][q])
            tm2["tidM" + str(i+1)].append(tm["tidM1"][q])
            taufrikm2["taufrikM" + str(i+1)].append(taufrikm["taufrikM" + str(i+1)][q])
            q += 1
        
        # print(q)
        
        while (q < len(tm["tidM"+str(i+1)])):
            om3["omegaM"+str(i+1)].append(om["omegaM"+str(i+1)][q])
            am3["alphaM" + str(i+1)].append(am["alphaM" + str(i+1)][q])
            tm3["tidM" + str(i+1)].append(tm["tidM" + str(i+1)][q])
            taufrikm3["taufrikM" + str(i+1)].append(taufrikm["taufrikM" + str(i+1)][q])
            q += 1

        
        """Plots"""
        #Hele turen
        hastighedsplot(tm["tidM" + str(i+1)],om["omegaM" + str(i+1)], 1)
        accelerationsplot(tm["tidM" + str(i+1)],am["alphaM" + str(i+1)], 2)
        friktionsmomentsplot(om["omegaM" + str(i+1)], taufrikm["taufrikM" + str(i+1)], 3)
        
        #Sammenligning med model
        hastighedsplot(tm["tidM" + str(i+1)],om["omegaM" + str(i+1)], 4)
        accelerationsplot(tm["tidM" + str(i+1)],am["alphaM" + str(i+1)], 5)

        #Kun før
        hastighedsplot(tm2["tidM" + str(i+1)],om2["omegaM" + str(i+1)], 6)
        accelerationsplot(tm2["tidM" + str(i+1)],am2["alphaM" + str(i+1)], 7)
        friktionsmomentsplot(om2["omegaM" + str(i+1)], taufrikm2["taufrikM" + str(i+1)], 8)

    """ Lineær regression og Printe tendenslinje """
    merge_omega_list_før += om2["omegaM"+str(i+1)]
    merge_omega_list_efter += om3["omegaM"+str(i+1)]
    merge_taufrik_list_før += taufrikm2["taufrikM" + str(i+1)]
    merge_taufrik_list_efter += taufrikm3["taufrikM" + str(i+1)]
    popt1, pcov1 = curve_fit(func,merge_omega_list_før,merge_taufrik_list_før)
    popt2, pcov2 = curve_fit(func,merge_omega_list_efter,merge_taufrik_list_efter)
    x1 = np.linspace(5,60,100)
    y1 = popt1[0]*x1+popt1[1]
    # x2 = np.linspace(5,60,100)
    # y2 = popt2[0]*x2+popt2[1]
    plt.plot(x1, y1, '-r', label='τ = ' + str(round(popt1[0],5)) + ' * ω + ' + str(round(popt1[1],5)))
    # plt.plot(x2, y2, '-r', label='τ = ' + str(round(popt2[0],5)) + ' * ω + ' + str(round(popt2[1],5)))

    
    #r^2-værdi
    y_pred_list = []
    for i in range(len(merge_omega_list_før)):
        y_pred_list.append(popt1[0]*merge_omega_list_før[i]+popt1[1])
    global r2
    r2 = r2_score(merge_taufrik_list_før,y_pred_list) ** (1/2)
    
    print('τ_friktion_før = ' + str(popt1[0]) + ' · ω + ' + str(popt1[1]))
    print('r^2-værdi = ' + str(r2**2))


    global slope_før
    slope_før = popt1[0]
    global intercept_før
    intercept_før = popt1[1]
    global slope_efter
    slope_efter = popt2[0]
    global intercept_efter
    intercept_efter = popt2[1]
    
       
    




def hjulberegner (hArr):
    for i in range(len(hArr)):
        """Vinkelhastighed- og acceleration af svinghjul"""
        oh["omegaS" + str(i+1)] = []
        oh2["omegaS" + str(i+1)] = []
        oh3["omegaS" + str(i+1)] = []
        ah["alphaS" + str(i+1)] = []
        ah2["alphaS" + str(i+1)] = []
        ah3["alphaS" + str(i+1)] = []
        tausys["tausys" + str(i+1)] = []
        svinginerti["svinginerti" + str(i+1)] = []

        th["tidS" + str(i+1)] = hArr[i]['time after start [s]'].tolist()
        th2["tidS" + str(i+1)] = []
        th3["tidS" + str(i+1)] = []
        rh["rotnumS" + str(i+1)] = hArr[i]['rotation number'].tolist()

        omegaberegner(th["tidS" + str(i+1)], oh["omegaS"+str(i+1)])
        alphaberegner(th["tidS" + str(i+1)], oh["omegaS"+str(i+1)], ah["alphaS"+str(i+1)])

        """Afskæring af data"""
        j = len(ah["alphaS" + str(i+1)])-2
        while ah["alphaS" + str(i+1)][j+1] < 0:
            j -= 1
        global tilstandsskifteS
        tilstandsskifteS = j

        q = 0
        while (q+1) < j:

            oh2["omegaS"+str(i+1)].append(oh["omegaS"+str(i+1)][q])
            ah2["alphaS" + str(i+1)].append(ah["alphaS" + str(i+1)][q])
            th2["tidS" + str(i+1)].append(th["tidS1"][q])
            q += 1
        q += 3
        while (q < len(oh["omegaS"+str(i+1)])):
            oh3["omegaS"+str(i+1)].append(oh["omegaS"+str(i+1)][q])
            ah3["alphaS" + str(i+1)].append(ah["alphaS" + str(i+1)][q])
            th3["tidS" + str(i+1)].append(th["tidS" + str(i+1)][q])
            q += 1

        """ Systemets moment og svinghjulets inertimoment"""
        for j in range(len(oh2["omegaS" + str(i+1)])):
            tausys["tausys" + str(i+1)].append(tau_flaske - (slope_før * oh2["omegaS" + str(i+1)][j] + intercept_før))
        # for j in range(len(oh3["omegaS" + str(i+1)])):
        #     tausys["tausys" + str(i+1)].append(slope_efter * oh3["omegaS" + str(i+1)][j] + intercept_efter)
        for j in range(len(ah2["alphaS" + str(i+1)])):
            svinginerti["svinginerti" + str(i+1)].append(tausys["tausys" + str(i+1)][j]/ah2["alphaS" + str(i+1)][j])
        # for j in range(len(ah3["alphaS" + str(i+1)])):
        #     svinginerti["svinginerti" + str(i+1)].append(tausys["tausys" + str(i+1)][j]/ah3["alphaS" + str(i+1)][j])

        #Hele turen
        hastighedsplot(th["tidS" + str(i+1)],oh["omegaS" + str(i+1)], 9)
        accelerationsplot(th["tidS" + str(i+1)],ah["alphaS" + str(i+1)], 10)

        #Kun før
        hastighedsplot(th2["tidS" + str(i+1)],oh2["omegaS" + str(i+1)], 11)
        accelerationsplot(th2["tidS" + str(i+1)],ah2["alphaS" + str(i+1)], 12)
        svinginertiplot(th2["tidS" + str(i+1)], svinginerti["svinginerti" + str(i+1)], 13)
        global MaaltInerti
        MaaltInerti = (np.mean(svinginerti["svinginerti" + str(i+1)]))
    print("Inertimoment af svinghjulet bliver " + str(np.mean(svinginerti["svinginerti" + str(i+1)])))
    
    


    
    
    

        
        
        


def omegaberegner(tid, omega):
    for i in range(len(tid)):
        if i == 0:
            omega.append((2 * math.pi)/tid[i])
        else:
            omega.append((2 * math.pi)/(tid[i]-tid[i-1]))


def alphaberegner(tid, omega, alpha):
    for i in range(len(omega)):
        if i == 0 :
            alpha.append(omega[i]/(tid[i] * 2))
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
    plt.legend(["Omega1","Omega2","Omega3","Model"])
    plt.xlabel('Tid [s]')
    plt.ylabel('ω [rad/s]')
    # plt.show

def accelerationsplot (tid, alpha, k):
    plt.figure(k).suptitle("Vinkelacceleration over for tid")
    plt.plot(tid,alpha)
    plt.legend(["Alpha1","Alpha2","Alpha3","Model"])
    plt.xlabel('Tid [s]')
    plt.ylabel('α [rad/s²]')
    # plt.show

def friktionsmomentsplot (vinkelhastighed, friktionsmoment, k):
    plt.figure(k).suptitle("Friktionsmoment over for vinkelhastighed")
    plt.plot(vinkelhastighed,friktionsmoment)
    plt.legend(["τ_frik1","τ_frik2","τ_frik3","Regression"])
    plt.xlabel('ω [rad/s]')
    plt.ylabel('τ_frik [N*m]')
    plt.show

def svinginertiplot (tid, svinginerti, k):
    plt.figure(k).suptitle("Inertimoment over for tid")
    plt.plot(tid,svinginerti)
    plt.legend(["τ_1","τ_2","τ_3","τ_4"])
    plt.xlabel('t [s] ')
    plt.ylabel('I [kg · s²]')
    # plt.show



murstensberegner(murstenArr)

hjulberegner(hjulArr)

MurstensModel()
SvinghjulsModel()
