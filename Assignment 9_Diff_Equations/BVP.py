import matplotlib.pyplot as plt
import my_lib as mb
import math



def dvdx(y,v,x):
    return 2*y
def dydx(y,v,x):
    return v


dt = 0.05
f = 1

l = -2
Y,V,X = mb.RK4_coup(dydx, dvdx, 1.2,l,0,dt,f)
yl = Y[19]
print(yl)
plt.plot(X,Y,'.',markersize=3)

h = -1
Y,V,X = mb.RK4_coup(dydx, dvdx,1.2,-1,0,dt,f)
yh = Y[19]
print(yh)
plt.plot(X,Y,'.',markersize=3)

g = mb.interpol(0.9,l,yl,h,yh)
print(g)
y,v,x = mb.RK4_coup(dydx, dvdx, 1.2, g,0,dt,f)
print(x[19])

plt.plot(x,y,'.',markersize=3)
# plt.plot()


G = [ele*0.005 for ele in range(0,200)]
H = [(0.157*math.exp((2**0.5)*ele)) + (1.043*math.exp(-1*(2**0.5)*ele)) for ele in G ]

plt.plot(G,H)
plt.show()
