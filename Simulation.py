#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 09:28:16 2020

@author: Tais
"""

import matplotlib.pyplot as plt
import Samlet

#Startvariabler
Start_hastighed = 700 #m/s
Start_omega = Start_hastighed / Samlet.r_sving
tid_list = [0]
omega_list = [Start_omega]
Energigraense = 100000 # Joule
alpha_list = [(-Samlet.slope_efter*Start_omega)/Samlet.MaaltInerti]
tabtenergi_list = [0]
t_interval = 1/10
energi_list = []

i = 0

while 1/2 * Samlet.MaaltInerti * omega_list[-1]**2 > Energigraense:
    tid_list.append(i*t_interval)
    omega_list.append(omega_list[i-1] + alpha_list[i-1] * (t_interval))
    alpha_list.append((-(Samlet.slope_efter * omega_list[i] + Samlet.slope_efter))/Samlet.MaaltInerti)
    i += 1

j = 0

for j in range(len(omega_list)):
    energi_list.append(1/2 * Samlet.MaaltInerti * omega_list[j]**2)
    
k = 0

for k in range(1,len(energi_list)):
    tabtenergi = energi_list[k-1] - energi_list[k] + tabtenergi_list[k-1]
    tabtenergi_list.append(tabtenergi)

#Plot
fig, ax1 = plt.subplots()
ax1.set_xlabel('Tid [s]')
ax1.set_ylabel('Kinetisk Energi [J]')
ax1.plot(tid_list,energi_list, color="blue")
ax1.tick_params(axis="y")
ax1.set_title('Energi over tid') 

fig, ax1 = plt.subplots()
ax1.set_xlabel('ω [rad/s]')
ax1.set_ylabel('Kinetisk Energi [J]')
ax1.plot(omega_list,energi_list, color="red")
ax1.tick_params(axis="y")
ax1.set_title('Kinetisk Energi over Vinkelhastighed')

fig, ax1 = plt.subplots()
ax1.set_xlabel('Tid [s]')
ax1.set_ylabel('ω [rad/s]')
ax1.plot(tid_list,omega_list, color="red")
ax1.tick_params(axis="y")
ax1.set_title('Vinkelhastighed over tid')

fig, ax1 = plt.subplots()
ax1.set_xlabel('Tid [s]')
ax1.set_ylabel('α [rad/s²]')
ax1.plot(tid_list,alpha_list, color="red")
ax1.tick_params(axis="y")
ax1.set_title('Vinkelacceleration over tid')

print("Det tager" , tid_list[-1], "sekunder før den kinetiske energi i hjulet"\
      " er under" , Energigraense, "Joule når hjulets hastighed ved kanten er", Start_hastighed, "m/s")
print("Systemet skal tilføres", tabtenergi_list[int(1/t_interval)], "W for at holde på den kinetiske energi")