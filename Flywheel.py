import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import stats
from sklearn.linear_model import LinearRegression

#Data Initialisers

data_sving = pd.read_csv("s11_1402.csv", sep = ",")

rotnum = data_sving['rotation number'].tolist()
tid = data_sving['time after start [s]'].tolist()



#Stang
m_stang = 0.2244 #kg
r_stang = 0.005 #m
I_stang = 1/2 * m_stang * r_stang ** 2 #kg * m^2

#Trisse
m_trisse = 0.015 #kg
r_trisse = 0.01 #m
I_trisse = m_trisse * r_trisse ** 2 #kg * m^2

#Plexi
N_plexi = 3
m_plexi = 0.057 #kg
r_plexi = 0.072 #m
I_plexi = N_plexi * m_plexi * r_plexi ** 2 #kg * m^2


#Møtrik
m_møtrik = 0.01 #kg
r1_møtrik = 0.005 #m
r2_møtrik = 0.008425 #m
N_møtrik = 11
I_møtrik = 1/2 * m_møtrik * (r1_møtrik ** 2 + r2_møtrik ** 2) * N_møtrik #kg * m^2

#Spændeskive
m_skive = 0.002 #kg
r1_skive = 0.005 #m
r2_skive = 0.009875 #m
N_skive = 9
I_skive = 1/2 * m_skive * (r1_skive ** 2 + r2_skive ** 2) * N_skive #kg * m^2

#Flaske
m_flaske = 2.078 #kg
g = 9.82 #N/kg
F_flaske = m_flaske * g #N
r_trisse = 0.01 #m
r_akse = 0.005 #m

#Svinghjul
m_sving = 4.980 #kg
r_sving = 0.4 #m

#Vinkelhastighed
omega_list = [] #m/s
for i in range(len(tid)):
    if i == 0:
        omega = (2 * math.pi)/(tid[i])
    else:
        omega = (2 * math.pi)/(tid[i]-tid[i-1])

    omega_list.append(omega)

#Vinkelacceleration
alpha_list = [] #m/s^2
for i in range(len(omega_list)):
    if i == 0:
        alpha = omega_list[i]/tid[i]
    else:
        alpha = (omega_list[i]-omega_list[i-1])/(tid[i]-tid[i-1])

    alpha_list.append(alpha)

omega_list2 = []
alpha_list2 = []
tid2 = []
i = 0
while alpha_list[i+1] > 0:
    omega_list2.append(omega_list[i])
    alpha_list2.append(alpha_list[i])
    tid2.append(tid[i])
    i += 1

#Systemets masse
m_system = m_sving + m_stang + m_trisse + m_plexi * N_plexi\
+ m_møtrik * N_møtrik + m_skive * N_skive + m_flaske #kg


#Normalkraft
N = m_system * g

#Friktionskoefficient
mu_list = []
a = 0.01206131003562931
b = 0.3997764794714707

for i in range(len(omega_list)):
    mu_list.append(a * omega_list[i] + b)


#Friktionskraft
frikkraft_list = []
for i in range(len(mu_list)):
    frikkraft_list.append(mu_list[i] * N)

#Friktionsmoment
taufrik_list = []
for i in range(len(frikkraft_list)):
    taufrik_list.append(frikkraft_list[i] * r_akse)





#Moment af flaske
tauflaske_list = []
for i in range(len(frikkraft_list)):
    tauflaske_list.append(F_flaske * r_trisse)


#Moment af Svinghjul
tausving_list = []
for i in range(len(tauflaske_list)):
    tausving_list.append(tauflaske_list[i] - taufrik_list[i])

#Inertimoment af svinghjul
svinginerti_list = []
for i in range(len(alpha_list2)):
    svinginerti_list.append(tausving_list[i]/alpha_list2[i])


'''Plots - før flaske'''

#Vinkelhastighed over for tid
fig2, ax2 = plt.subplots()
ax2.set_xlabel('t [s]')
ax2.set_ylabel('ω [rad/s]')
ax2.plot(tid,omega_list, color="blue")
ax2.tick_params(axis="y")
ax2.set_title('Vinkelhastighed over for tid - Svinghjul Før')

#Vinkelacceleration over for tid
fig, ax1 = plt.subplots()
ax1.set_xlabel('t [s]')
ax1.set_ylabel('α [rad/s^2]')
ax1.plot(tid2,alpha_list2, color="red")
ax1.tick_params(axis="y")
ax1.set_title('Vinkelacceleration over for tid - Svinghjul Før')

#Inertimoment over for tid
fig3, ax3 = plt.subplots()
ax3.set_xlabel('t [s]')
ax3.set_ylabel('I [kg * m^2]')
ax3.plot(tid2,svinginerti_list, color="silver")
plt.ylim(-.5, 1)
ax3.set_title('Inertimoment over for tid - Svinghjul Før')

#Beregnet Inertimoment
x = np.linspace(5,70,100)
y = 0*x + 0.5817162071
plt.plot(x, y, '-r', label='Beregnet inertimoment')

plt.legend(loc='upper left')


plt.show()
