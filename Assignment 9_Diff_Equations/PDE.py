import my_lib as mb
import math



lx = 2
lt = 4
nx = 20
nt = 5000
hx = lx/nx
ht = lt/nt

import matplotlib.pyplot as plt

def PDE(X,profile,nx,nt,hx,ht):
    a = ht/(hx*hx)

    for j in range(0,nt+1):
        if j <= 5000:
            if j == 0 or j == 10 or j == 20 or j == 50 or j == 100 or j == 200 or j == 500 or j == 1000:
                plt.plot(X,profile,label="t = "+str(j*ht))

        for k in range(0,nx+1): # matrix multiplication A*Vj
            if k==1:
                profile[k] = profile[k+1]*a + profile[k]*(1-2*a)
            if k==nx:
                profile[k] = profile[k-1]*a + profile[k]*(1-2*a)
            elif k!=nx and k!=0:
                profile[k] = profile[k-1]*a + profile[k]*(1-2*a) + profile[k+1]*a
    plt.legend()
    plt.show()

PDE(a,nx,nt,hx,ht)
