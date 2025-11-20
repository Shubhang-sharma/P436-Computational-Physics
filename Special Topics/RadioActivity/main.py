seed = 10
import math
import matplotlib.pyplot as plt
# In a time of 60 minutes how many times you are seeing 1 decay, .....

def myrand(a=1103515245,c = 12345,m = 32768): # returns random numbers in range(0,1)
    global seed
    seed = ((a*seed + c)%m)/m
    return seed

Ta,Tb,dt,T,na,nb,nc = 20,30,2,0,500,0,0

la = math.log(2)/Ta
lb = math.log(2)/Tb

Na,Nb,Nc,TR = [500],[0],[0],[0]

while T <= 10*Ta:
    T+=dt
    for i in range(na):
        if myrand() <= la*dt:
            na-=1
            nb+=1
    Na.append(na)
    Nb.append(nb)
    TR.append(T)
plt.plot(TR,Na,".")
plt.plot(TR,Nb,".")


while T <= 10*Ta:
    T+=dt
    for i in range(na):
        if myrand() <= la*dt:
            na-=1
            nb+=1
    for i in range(nb):
        if myrand() <= lb*dt:
            nb-=1
            nc+=1
    Na.append(na)
    Nb.append(nb)
    Nc.append(nc)
    TR.append(T)

plt.plot(TR,Na,".",label="Na")
plt.plot(TR,Nb,".",label="Nb")
plt.plot(TR,Nc,".",label="Nc")
plt.legend()
plt.plot()





plt.show()
