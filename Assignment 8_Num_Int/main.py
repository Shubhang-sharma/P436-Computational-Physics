import my_lib as mb
import math
import matplotlib.pyplot as plt
output = open('output.txt','w')
Res = []

#Question-1:
Res.append("Question-1 \n")

# ----------
def f(x):
    return (1+1/x)**0.5
Res.append("      Midpoint  Trapezoidal   Simpson \n")
for N in range(10,31,10):
    Res.append("N = "+str(N)+"  "+str(round(mb.midpoint(f,1,4,N),4)) + "       ")
    Res.append(str(round(mb.trapezoidal(f,1,4,N),4)) + "       ")
    Res.append(str(round(mb.simpson(f,1,4,N),4)) + " \n")
#-----------


#Question-2:
Res.append("\n Question-2 \n")

#----------
def g(x):
    return math.sin(x)**2

n = 5000
int,sg,FN,SG,N = mb.monte(g,-1,1,n)
print(sg)
Res.append('Integral via Monte Carlo Method= ' + str(round(int,4)))

plt.plot(N,FN,'.',markersize=1.5,label = "FN")
plt.plot(N,SG,'.',markersize=1.5,label = "$\sigma_{f}$")
plt.plot(N,[int for i in range(len(N))],'.',markersize=0.8,label = "FN =" +str(round(int,4))+" at N = " +str(n))
plt.xlabel('N')
plt.title("Monte Carlo Sampling with #Random Numbers (N) = 10 to "+ str(n)+"\n \n Numerical Integration value = " + str(round(int,4)))
plt.legend(loc = "upper right")
plt.show()
#------------------


#Question-3:
Res.append("\n Question-3 \n")

#-----------
def lin_density(x):
    return x*x
def moment1(x):
    return lin_density(x)*x

com = mb.simpson(moment1,0,2,10)/mb.simpson(lin_density,0,2,10)

Res.append('Centre Of Mass = ' + str(round(com,5)) + ' units from the reference end')
#-----------


#OUTPUT------------------------------------
for i in range(len(Res)):
    # for j in range(len(Res[i])):
    output.writelines(Res[i])
output.close()
