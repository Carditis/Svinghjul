
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import stats


#Data Initialisers

data_sten = pd.read_csv("m11_1423.csv", sep = ",")


rotnum = data_sten['rotation number'].tolist()
tid = data_sten['time after start [s]'].tolist()




#Fjerner data for når flasken er sluppet

for i in range(7):
    del rotnum[-1]
    del tid[-1]




#Intertimoment 
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

I = I_sten + I_stang + I_trisse + I_plexi + I_møtrik + I_skive #kg * m^2


#Svinghjul
m_svinghjul = 4.980 #kg

#Vinkelhastighed
omega_list = [] #m/s
for i in range(len(tid)):
    if i == 0:
        omega = (2*math.pi)/tid[i]
    else:
        omega = (2 * math.pi)/(tid[i]-tid[i-1])
        
    omega_list.append(omega)
    
#Vinkelacceleration
alpha_list = [] #m/s^2
for i in range(len(omega_list)):
    if i == 0:
        alpha = (omega_list[i])/tid[i]
    else:
        alpha = (omega_list[i]-omega_list[i-1])/(tid[i]-tid[i-1])

    alpha_list.append(alpha)
        
   
#Moment 
tausystem_list = [] #N*m/kg
for i in range(len(alpha_list)):
    tau = I * alpha_list[i]
    
    tausystem_list.append(tau)
    
    
#Flaske
m_flaske = 2.078 #kg
g = 9.82 #N/kg
F_flaske = m_flaske * g #N
r_trisse = 0.01 #m
r_akse = 0.005 #m
tau_flaske = F_flaske * r_trisse #N*m/kg


#Friktion
taufrik_list = [] #N*m/kg
for i in range(len(tausystem_list)):
    taufrik_list.append(-tausystem_list[i]+tau_flaske)

#tau_flaske + tau_frik = sum(tau_list)
    
    
#Friktionskraft
frikkraft_list = [] #N
for i in range(len(taufrik_list)):
    frikkraft_list.append(taufrik_list[i]/r_akse)
  
    
    
#Normalkraft
    m_system = m_sten + m_stang + m_trisse + m_plexi * N_plexi\
    + m_møtrik * N_møtrik + m_skive * N_skive + m_flaske #kg
    N = m_system * g #N
    
    
#Friktionskoefficient
mu_list = [] 
for i in range(len(frikkraft_list)):
    mu_list.append(frikkraft_list[i]/N)
    



""" Vinkelhastighed over for tid """
            
fig3, ax3 = plt.subplots()
ax3.set_xlabel('Tid [t]')
ax3.set_ylabel('Vinkelhastighed [ω]')
ax3.plot(tid, omega_list, color="blue")
ax3.tick_params(axis='y')
ax3.set_title('Vinkelhastighed over for Tid')



#Vinkelacceleration over for tid
fig2, ax2 = plt.subplots()
ax2.set_xlabel('Tid [t]')
ax2.set_ylabel('Vinkelacceleration [α]')
ax2.plot(tid,alpha_list, color="blue")
ax2.tick_params(axis="y")
ax2.set_title('Vinkelacceleration over for Tid')


'''Systemets moment over for tid'''
fig5, ax5 = plt.subplots()
# Dernæst beskrives "ax1"'s x og y akser
ax5.set_xlabel('Tid [t]')
ax5.set_ylabel('Moment')
# Her plottes (V,I) data som punkter
ax5.scatter(tid, tausystem_list, color="blue")
# Her farves teksten der beskrvier y-aksen
ax5.tick_params(axis='y')
ax5.set_title('Systemets Moment over for Tid')

#Friktionskoefficient over for vinkelhastighed
fig, ax1 = plt.subplots()
# Dernæst beskrives "ax1"'s x og y akser
# Her plottes (V,I) data som punkter
ax1.plot(omega_list, mu_list, color="blue")
ax1.set_title('Friktionskoefficient over for vinkelhastighed')
x = np.linspace(5,40,100)
y = 0.01206131003562931*x+0.3997764794714707         
plt.plot(x, y, '-r', label='y=0.012x+0.400\
         r^2 = 0.67')
plt.xlabel('Vinkelhastighed', color='#1C2833')
plt.ylabel('Friktionskoefficient', color='#1C2833')
plt.legend(loc='upper left')  

# Her farves teksten der beskrvier y-aksen
ax1.tick_params(axis='y')

plt.show()


'''Lineær Regression på friktionskoefficient'''
slope, intercept, r_value, p_value, std_err = stats.linregress(omega_list,mu_list)

a = slope
b = intercept



print("y = " + str(slope) + " * x + " + str(intercept))
print("r^2 = " + str(r_value**2))






#plt.axis((15,x2,y1,y2))











    
    








