import random

f = open("weight.csv","w")

N=2000 #The dimension of your directionaless graph

for i in range(2,N+1):
    j=1
    while j<i:
        w_ij= random.choice((-1,1))
        f.write(str(i) + "," + str(j) +"," + str(w_ij) + "\n")
        j += 1

f.close