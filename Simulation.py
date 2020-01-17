#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 09:28:16 2020

@author: Tais
"""

import matplotlib.pyplot as plt

#Startvariabler
Start_hastighed = 700 #m/s
r_sving = 0.4
Start_omega = 700 * 0.4
tid_list = [0]
omega_list = [Start_omega]
slope = 0.002874918531194114
intercept = 0.09479051847376596
I = 0.6483255722121865
Energigraense = 20000
alpha_list = [(-intercept)/I]

i = 0

while 1/2 * I * omega_list[-1]**2 > Energigraense:
    tid_list.append(i/10)
    omega_list.append(omega_list[i-1] + alpha_list[i-1] * (1/10))
    alpha_list.append((-(slope * omega_list[i] + intercept))/I)
    i += 1

energi_list = []

j = 0

for j in range(len(omega_list)):
    energi_list.append(1/2 * I * omega_list[j]**2)

#Plot
fig, ax1 = plt.subplots()
ax1.set_xlabel('Tid [s]')
ax1.set_ylabel('Kinetisk Energi [J]')
ax1.plot(tid_list,energi_list, color="blue")
ax1.tick_params(axis="y")
ax1.set_title('Energi over tid') 

# fig1, ax1 = plt.subplots()
# ax1.set_xlabel('Vinkelhastighed [omega]')
# ax1.set_ylabel('Kinetisk Energi [J]')
# ax1.plot(omega_list,energi_list, color="red")
# ax1.tick_params(axis="y")
# ax1.set_title('Kinetisk Energi over Vinkelhastighed')
print("Det tager " , tid_list[-1], " sekunder før den kinetiske energi i hjulet "\
      " er under " , Energigraense, " Joule når hjulets hastighed ved kanten er", Start_hastighed, " m/s")