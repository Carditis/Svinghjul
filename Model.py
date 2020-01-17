import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import stats
from sklearn.linear_model import LinearRegression
from Mursten import slope, intercept

'''Inertimomenter'''
#Sten
a = 0.088 #m
b = 0.224 #m
m_sten = 2.224 #kg
I_sten = m_sten * 1/12 * (a*a + b*b) #kg * m^2

#Stang
m_stang = 0.2244 #kg
r_stang = 0.005 #m
I_stang = 1/2 * m_stang * r_stang * r_stang #kg * m^2

#Trisse
m_trisse = 0.015 #kg
r_trisse = 0.01 #m
I_trisse = m_trisse * r_trisse * r_trisse #kg * m^2
    
#Plexi
N_plexi = 3
m_plexi = 0.057 #kg
r_plexi = 0.072 #m
I_plexi = N_plexi * m_plexi * r_plexi * r_plexi #kg * m^2


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



'''Momenter'''
m_flaske = 2.078 #kg
g = 9.82 #N/kg
r_trisse = 0.01 #m
r_akse = 0.005 #m

tau_flaske = m_flaske * g * r_trisse

'''Lister'''
omega_list = []
omega_list.append(0)
alpha_list = []
alpha_list.append((tau_flaske - intercept)/I)
tid_list = []
tid_list.append(0)

L_snor = 7.2
i = 1

while 2 * math.pi * r_trisse * i < L_snor:
    tid_list.append(i/10)
    omega_list.append(omega_list[i-1] + alpha_list[i-1] * (1/10))
    alpha_list.append((tau_flaske - (slope * omega_list[i] + intercept))/I)  
    i += 1

while omega_list[i-1] > 0:
    tid_list.append(i/10)
    omega_list.append(omega_list[i-1] + alpha_list[i-1] * (1/10))
    alpha_list.append((-(slope * omega_list[i] + intercept))/I)
    i += 1
         
         
'''Plots'''
#Vinkelhastighed over for tid
fig, ax1 = plt.subplots()
ax1.set_xlabel('Tid [s]')
ax1.set_ylabel('Vinkelhastighed [rad/s]')
ax1.plot(tid_list,omega_list, color="blue")
ax1.tick_params(axis="y")
ax1.set_title('Vinkelhastighed over for tid - Model Før/Efter') 



#Vinkelacceleration over for tid
fig2, ax2 = plt.subplots()
ax2.set_xlabel('Tid [s]')
ax2.set_ylabel('Vinkelacceleration [rad/s^2]')
ax2.plot(tid_list,alpha_list, color="red")
ax2.tick_params(axis="y")
ax2.set_title('Vinkelacceleration over for tid - Model Før/Efter')

   

