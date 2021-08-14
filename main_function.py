import random
import numpy as np
import pandas as pd
import math 
from random import choice
import math
import os.path
from time import process_time
import threading

#number of verties
N=2000

#timer
total_time = 0
time_construct = 0

#creative a random spin array
spin_list = []
for i in range(N):
    spin_list.append(random.choice((-1,1)))
spin_array_cons = np.array(spin_list)

def IsingEnergy(spin_array, weight):
    energy = 0
    for i in range(len(weight)):
        particle_1 = weight[i][0]
        particle_2 = weight[i][1]
        w_12 = weight[i][2]

        energy += w_12*spin_array[particle_1 - 1]*spin_array[particle_2 - 1]
    return(energy)

def partial_changing_isingenergy(spin_array, weight, random_int):
    a = random_int
    delta_energy=0
    if(a == 0):
        for i in range(N-1):
            t = int(i*i/2)
            w_12 = weight[t][2]
            particle_1 = weight[t][0]
            particle_2 = weight[t][1]
            delta_energy += w_12 * spin_array[particle_1 - 1] * spin_array[particle_2 - 1]
        return (delta_energy)

    else:
        delta_energy_1 = 0
        delta_energy_2 = 0
        k_1 = int((a*a-a)/2 + a - 1)
        w = weight[k_1][2]
        p_1 = weight[k_1][0]
        p_2 = weight[k_1][1]
        delta_energy_1 = w * spin_array[p_1 - 1] * spin_array[p_2 - 1]
        if(a != N-1):
            for i in range(a+1,N):
                k_2 = int((i*i-i)/2 + a -1)
                w_12 = weight[k_2][2]
                particle_1 = weight[k_2][0]
                particle_2 = weight[k_2][1]
                delta_energy_2 += w_12 * spin_array[particle_1 -1] * spin_array[particle_2 - 1]
        delta_energy = delta_energy_1 + delta_energy_2
        return(delta_energy)

def SummationWeight(weight):
    sum_weight = 0
    for i in range(len(weight)):
        sum_weight += weight[i][2]
    return(sum_weight)

def num_cut(sum_weight,energy):
    cut = (sum_weight - energy)/2
    return(cut)

def flipping(spin_array,random_int,energy_0,weight,beta):
    spin_array[random_int] = -spin_array[random_int]
    Delta_energy = 2 * partial_changing_isingenergy(spin_array,weight,random_int)
    if(Delta_energy < 0 or math.exp(-beta * Delta_energy) >= np.random.uniform(0,1) ):
        energy_0 = energy_0 + Delta_energy
        accept_not = 1
    else:
        spin_array[random_int] = -spin_array[random_int]
        accept_not = 0
    return(energy_0,spin_array,accept_not)

def solve_system(spin_array,weight,beta_0,T,sum_weight,best_cut,energy_0,target_energy,time_limit):
    time_start = process_time()
    energy_0 = IsingEnergy(spin_array, weight)
    t=0
    while t < T+1: 
        random_int = random.randint(0,N-1)

        beta = beta_0*math.log(1+t/T)
        flip = flipping(spin_array,random_int,energy_0,weight,beta)
        energy_0 = flip[0]
        spin_array = flip[1]
        h = flip[2]
        cut = num_cut(sum_weight,energy_0)
        if(h == 1):
            t += 1
        if(target_energy != None):
            if(energy_0 <= target_energy):
                break
        if(time_limit != None):
            time_end_prime = process_time()
            delta_time_prime = time_end_prime - time_start
            if(delta_time_prime >= time_limit):
                break
        if(cut>best_cut):
            best_cut = cut
    time_end = process_time()
    delta_time = time_end - time_start
    return (energy_0,best_cut,delta_time,spin_array)

#import weight from weight.scv
file_name="weight.csv"
BASE_DIR=os.path.dirname(os.path.abspath(__file__))
file=os.path.join(BASE_DIR,file_name)
print(file)
data = pd.read_csv(file,header=None)

#parameters & variables
weight = data.values.tolist()
sum_weight = SummationWeight(weight)

beta_0 = 4.0
Time_scaling_factor = 4200
target_energy = -60000 #int, float number or None
time_limit = None #int, float number or None
num_try = 10

total_Hamiltonain = 0
best_cut = 0
best_energy = 0
total_cut = 0

for j in range(num_try):
    result = 0
    spin_array = spin_array_cons
    result = solve_system(spin_array,weight,beta_0,Time_scaling_factor,sum_weight,0,best_energy,target_energy,time_limit)
    if(result[0]<best_energy):
        best_energy = result[0]
    if(result[1]>best_cut):
        best_cut = result[1]
    total_Hamiltonain += result[0]
    total_cut += result[1]
    total_time += result[2]
    print(result[0],result[1],result[2])

avg_Hamiltonian = total_Hamiltonain/num_try
avg_cut = total_cut/num_try
avg_time = total_time/num_try
print("------------------------------------")
print("CONDITIONS")
print("------------------------------------"+"\n")
print("Trials:",num_try)
print("Number of vertices: 2000")
print("Temperature scaling factor beta_0:",beta_0)
print("Time scaling factor T:",Time_scaling_factor)
print("Target energy:", target_energy)
print("Time limit:",time_limit,"s")
print("------------------------------------")
print("RESULT")
print("------------------------------------"+"\n")
print("Best cut:",best_cut)
print("Best Ising energy:",best_energy)
print("Average cut:",avg_cut)
print("Average processing time:",avg_time,"s")
print("Average energy finally reached:", avg_Hamiltonian)
print("------------------------------------")